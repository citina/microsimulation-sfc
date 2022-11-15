% This function will take in a state matrix and calculate the infection
% probabilities for each demographic group using the formula 
% P(infection) = (I/N) * S

function infection_probability = calc_infec_prob4(state_matrix, StateMatCols, mix_demog_group_def, mix_demog_col_struct, common_col_struct, artAdherence)

    % Find the people in the demographic group
    demog_mat_indices = find_demog_rows(state_matrix, StateMatCols, mix_demog_group_def, mix_demog_col_struct, common_col_struct);
    
    
    
    % If there are no people in the demographic group, return an infection
    % probability of 0
    if size(demog_mat_indices, 2) == 0
       infection_probability = 0;
       return;
    end

    % find all the transmitting people in the demographic group (I) -
    % defined by the following
        % patient is alive ()
        % status is not suceptible(not 0)
        % does not have treatment (treated is 0)
        % has treatment but does not adhere
    alive = find_indices(state_matrix, demog_mat_indices, StateMatCols.alive, '=', 1);

    % find total number of people in the demographic group (N) whoa re alive
    N = size(alive, 2);

    % if nobody alive in demographic group, add 0
    if N == 0
       infection_probability = 0;
       return;
    end

    has_disease = find_indices(state_matrix, alive, StateMatCols.status, '>', 0);
    has_disease_no_treatment = find_indices(state_matrix, has_disease, StateMatCols.treatment, '=', 0);
    has_disease_treatment = find_indices(state_matrix, has_disease, StateMatCols.treatment, '=', 1);
    has_disease_no_treatment_prep = find_indices(state_matrix, has_disease_no_treatment, StateMatCols.prep, '=', 1);
    has_disease_no_treatment_no_prep = find_indices(state_matrix, has_disease_no_treatment, StateMatCols.prep, '=', 0);

    I = size(has_disease_no_treatment_no_prep, 2) + size(has_disease_no_treatment_prep, 2) + size(has_disease_treatment,2)* (1-artAdherence);
    % I
    
    infection_probability = I/N;

    % infection_probability = 0;
    
end