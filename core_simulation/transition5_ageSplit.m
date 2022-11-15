% NEED TO DOCUMENT THE CHANGES I MADE!!!!!


function [updated_state_matrix, transitionTally] = transition5_ageSplit(transition,transition_path, state_matrix, StateMatCols, state_mat_demog_group_idx, future_state_param_names, calib_a1, calib_a2, calib_a3, calib_a4)
    
    % make a copy of the state matrix such that a tally can be made of the number of changes that occur
    state_matrix_copy = state_matrix;

    % Read in transition matrix
    [state_trans_mat, StateTransCols] = read_table(transition_path);

    
    % Find the people who are in the starting state within the
    % demographic group
    eligible_rows_demog = state_mat_demog_group_idx;


       
    % iterate through the rows of the transition matrix.
    % each row in the transition matrix has a different probability for
    % transition
    for initialState = 1:size(state_trans_mat,1)
        
        % reset possible eligible rows for each iteration of a new inital
        % state
        eligible_rows = eligible_rows_demog;

        % iterate through the field names in the transition
        % use this to find the final set of eligible rows
        for field = fieldnames(StateTransCols)'
            field_string = string(field);
            if(~contains(field_string,'state') && ~contains(field_string,'prob'))
                start_state_value = state_trans_mat(initialState, StateTransCols.(field_string));
                % check if fieldname contains age_min or age_max
                if(contains(field_string,'age_min'))
                    eligible_rows = find_indices(state_matrix, eligible_rows, StateMatCols.age, '>=', start_state_value);
                elseif (contains(field_string,'age_max'))
                    eligible_rows = find_indices(state_matrix, eligible_rows, StateMatCols.age, '<=', start_state_value);
                else
                    eligible_rows = find_indices(state_matrix, eligible_rows, StateMatCols.(field_string), '=', start_state_value); 
                end
            end
        end


        % Get all state level index
        state_idx = find(contains(fieldnames(StateTransCols), future_state_param_names(1)));

        % Get all prob index
        prob_idx = find(contains(fieldnames(StateTransCols), future_state_param_names(2)));

        
        % Look up probabilities from the table associated with this
        % particular initial state. 
        % Store in array to be used for comparisons
        % Look up associated states
        % store in array
        drawComp = ones(1, length(prob_idx));
        newState = ones(1, length(state_idx));
        
        for i = 1:length(prob_idx)
            drawComp(i) = state_trans_mat(initialState,prob_idx(i));
            newState(i) = state_trans_mat(initialState,state_idx(i));
        end

        


        % apply calibration constant to scale the probability applied here
        % calibration constant depends on the age of the individual based on a cutoff of 50. 
        % age index in state matrix is 1
        ageIndex = 1;
        % 50 is the place where we split age off of. Note that 50 is associated IS associated with our age strata
        ageCut_1 = 30;
        ageCut_2 = 50;
        ageCut_3 = 65;
        
        % each transition represents a single age year
        % we thus just need to check that the age for the elible rows is either old or young


        % only go through this step if have elgible rows
        if length(eligible_rows > 0)

            % this works because for death we are age specific
            eligibleAge = state_matrix_copy(eligible_rows(1), ageIndex);

            if eligibleAge < ageCut_1
                calibrationConstant = calib_a1;
            elseif eligibleAge < ageCut_2
                calibrationConstant = calib_a2;
            elseif eligibleAge < ageCut_3
                calibrationConstant = calib_a3;
            else
                calibrationConstant = calib_a4;
            end
        else
            % proceed as if no calibration constant if 0
            calibrationConstant = 1;
        end

        drawComp = drawComp *calibrationConstant;

        % draw comp cannot surpass 1 or the process does not work... this is why we run into issues with death
        drawComp = min(drawComp, 1);




        % reduce drawComp and newState to only entries were prob of transition is greater than 0
        keepIndex = drawComp > 0;
        drawComp = drawComp(keepIndex);
        newState = newState(keepIndex);


        % make the items remaining in drawComp effectively a cumulative freq
        for i = 2:length(drawComp)
            drawComp(i) = drawComp(:, i-1) + drawComp(:,i);
        end

        % take 1 minus the new cumulative freq and reverse the order
        % the state must also be flipped, but not 1 minus
        % this is because we want to transition if random draw is greater value 
        % in arraway associated with the state chagne. If greater, than the state will change
        % EX: states 1, 2, 3 with prob 0, .2, and .4
            % as a cumulative frequency vector, we will have [.2 .6] assicoated with states [2 3]
            % be reversing the order and taking the 1 minus of that, we have [.4 .8] associated with [3 2]

            % using a random draw and the iterative approach, there is a .4 chance no state change will happen. 
            % we then have a .8-.4 =.4 chance a change to state 3 will occur
            % after the transition to state 3 is checked, the transition to state 2 is checked
            % we have a 1-.8 = .2 probability of transition to state 2

            % note that because of the for loop, a value of .9 for the random draw would physically transition 
                % to state 3 but then again be transitioned to show state 2 as the end result of the loop. 
        
        drawComp = 1-fliplr(drawComp);
        newState = fliplr(newState);
      
        % create a random draw value for all eligible rows
        unif_prob = rand([1,length(eligible_rows)]);
        
        % iterate through each comparison of the cumulative frequencies to
        % determine 
        for i = 1:length(drawComp)
           
           
            to_transition = drawComp(i) <= unif_prob;
            
            % identify which transitions should occur
            state_mat_transition_rows = eligible_rows(to_transition);

            
            % transition people to appropriate stae
            state_matrix(state_mat_transition_rows, StateMatCols.(transition)) = newState(i);
        end
            
    end

    A = state_matrix_copy(:, StateMatCols.(transition)) ~= state_matrix(:, StateMatCols.(transition));


    transitionTally = sum(A);

    
    % if transition == "alive"
    %     transition_path
    %     marker = "NewDemogGroup"
    %     demogSize = length(eligible_rows_demog)
    %     drawComp
    %     transitionTally
    % end
    
    updated_state_matrix = state_matrix;
end
      
