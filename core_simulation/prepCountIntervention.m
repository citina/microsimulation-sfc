


function [updated_state_matrix, interventionTally] = prepCountIntervention(interventionPath, state_matrix, StateMatCols)


     % make a copy of the state matrix such that a tally can be made of the number of changes that occur
     state_matrix_copy = state_matrix;


     % Read in intervention plan
    [intervention_mat, interventionCols] = read_table(interventionPath);


    % get index for the intervetnion matrix (not same as state matrix)
    % this is hardcoded but we could use fieldnames (see transition.m)
    ageMin_index = 1;
    ageMax_index = 2;
    % we are not consider race in this intervention
    % race_index = 3;
    % status_index = 4;
    aware_index = 3;
    prep_index = 4;
    alive_index = 5;
    quant_index = 6;


    % iterate through the rows of intervention groups
    % each row specifies the number of people being impacted
    num_subgroups = size(intervention_mat,1);

    for i = 1:num_subgroups

        subgroup = intervention_mat(i,:);
        
        % parameters we are itnerested in
        ageMin = subgroup(ageMin_index);
        ageMax = subgroup(ageMax_index);
        % ignoring race and status index. status must be zero
        % race = subgroup(race_index);
        % status  = subgroup(status_index);
        aware = subgroup(aware_index);
        prep = subgroup(prep_index);
        alive = subgroup(alive_index);

        % number of prep
        numPrep = subgroup(quant_index);

        % get index of people intersted in
        aliveVec = state_matrix(:,StateMatCols.alive);
        ageVec = state_matrix(:,StateMatCols.age);
        raceVec = state_matrix(:,StateMatCols.race);
        statusVec = state_matrix(:,StateMatCols.status);
        awareVec = state_matrix(:,StateMatCols.aware);
        prepVec = state_matrix(:,StateMatCols.prep);

        % get the index of all individuals who fit criteria
        boolVec = (ageVec >= ageMin & ageVec <= ageMax & statusVec == 0 & awareVec == aware & prepVec == prep & aliveVec == alive);

        % get indices of interest
        indices = find(boolVec);

        length(indices);

        % if the number of people who can start PrEP is smaller than desired amount to give
            % set the num on prep to the min
        numPrep = min([numPrep, length(indices)]);

        % randomly sample the number of individuals without replacement
        givePrep = datasample(indices, numPrep, 'Replace', false);

        startPrep = ones(numPrep, 1);

        % change the indices of the people identified
        state_matrix(givePrep, StateMatCols.prep) = startPrep;

    end


    % double check that the proper number of changes were made

    A = state_matrix_copy(:, StateMatCols.prep) ~= state_matrix(:, StateMatCols.prep);

    interventionTally = sum(A);

    % display("This is the number of PrEP we added via intervention")
    % interventionTally

    updated_state_matrix = state_matrix;

end


