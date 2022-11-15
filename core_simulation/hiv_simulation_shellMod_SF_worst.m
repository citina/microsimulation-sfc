% tic % start timer

% paths
code_loc = pwd;
state_matrices_path = sim_dataDir;
addpath([fileparts(code_loc), '/global_fun']);
addpath([fileparts(code_loc), '/input']);

% input files (initial pop and workbook)
init_yr = 2012;
initialPopFile = sim_initialPop;
inputWorkbook = sim_inputFile;

% A flag to test policy. The value of the flag is the year the policy starts
policyFlag = interventionStart; % No interventions

% A flag for type 2 policy (0 or 1: intervention or not)
policyFlag_type2 = interventionFlag2;

% define input parameters and paths of inputs
hiv_simulation_parameters_SF;
prepAdherence_default = prepAdherence;
prepOff_transition_path_default = prepOff_transition_path;
treatmentOn_transition_path_default = treatmentOn_transition_path;
aware1_transition_path_default = aware1_transition_path;
aware2_transition_path_default = aware2_transition_path;
aware3_transition_path_default = aware3_transition_path;

% Read in the initial population
% Assumption: the file has all the people pre-populated with characteristics
[state_matrix, StateMatCols] = read_table(init_pop_file);

% Read in all possible demographic groups
% age min and max MUST be listed first in excel sheet
[demog_table, DemogTblCols] = create_demog_groups(demog_var_def_file);

% Read in demographic groups specified in the mixing matrix
% age min and max MUST be listed first in excel sheet
[mixing_table, MixingTblCols] = create_demog_groups(mixing_mat_def_file);

% Read in the mixing matrix
mixing_matrix = importdata(mixing_mat_file);

% read in number of partners vector
numPartners_default = importdata(num_partners_file);

% empty matrix to store tally
finalTransitionTally = [];

% clear current matrices in folder
cd(sim_dataDir);
delete *.mat

% save initial population as a state matrix
matrix_name = int2str(0);
save(matrix_name, 'state_matrix');

%% covid
% covid death
p_covid_d = zeros(1, T);
p_covid_d_plwh = zeros(1, T);
p_covid_d(yr_covid) = [0.000246, repelem(0.0003144, length(yr_covid)-1)];
p_covid_d_plwh(yr_covid) = [0.0004076, repelem(0.0005209, length(yr_covid)-1)];
p_covid_d_race = [0.094, 0.204, 0.289, 0.413]; %B, H, W, O
p_covid_d_age = [0, 0.004, 0.015, 0.038, 0.073, 0.145, 0.176, 0.549];

covid_age_bucket = [15, 19, 20, 29, ...
                    30, 39, 40, 49, ...
                    50, 59, 60, 69, ...
                    70, 79, 80, 99];
                
prepAdherence_worst = [0.024096386 0.33686747 0.63903614];

%% For each time period
cd(code_loc)
for t = 1:T
    % update prep adoption rate 
    if t == 1
       prepOn_transition_path = prepOn_transition_path_2013; 
    end
    if t == 2
       prepOn_transition_path = prepOn_transition_path_2014; 
    end
    if t == 3
       prepOn_transition_path = prepOn_transition_path_2015; 
    end
    if t == 4
       prepOn_transition_path = prepOn_transition_path_2016; 
    end
    if t == 5
       prepOn_transition_path = prepOn_transition_path_2017; 
    end 
    if t >= 6 % yr 2018
       prepOn_transition_path = prepOn_transition_path_2018;
    end

    %% covid prep -start-
    if ismember(t, yr_prep_init) && prep_init_flg
       prepOn_transition_path = file_common_part + 'input/transitions/prepOn_transition_worst.csv';
    end
    if ismember(t, yr_prep_dis) && prep_dis_flg
       prepOff_transition_path = file_common_part + 'input/transitions/prepOff_transition_worst.csv';
       prepAdherence = prepAdherence_worst;
    else
       prepAdherence = prepAdherence_default; %reset it back to non covid prep adherence
       prepOff_transition_path = prepOff_transition_path_default;
    end   
    %% covid prep -end-
    %% covid treatment -start-
    if ismember(t, yr_vls) && vls_flg
       treatmentOn_transition_path = file_common_part + 'input/transitions/treatmentOn_transition_worst.csv';
    else
       treatmentOn_transition_path = treatmentOn_transition_path_default;
    end   
    %% covid treatment -end- 
    %% covid num of partners -start-
    partners_re = 0.125;
    if ismember(t, yr_par) && par_flg 
        numPartners = numPartners_default * (1-partners_re);
    else
        numPartners = numPartners_default;
    end
    %% covid num of partners -end-
    
    %% covid death setting -start-  
    alive = state_matrix(:,StateMatCols.alive);
    age = state_matrix(:,StateMatCols.age);
    race = state_matrix(:,StateMatCols.race);
    status = state_matrix(:,StateMatCols.status);

    % covid death assignment
    if p_covid_d(t) ~= 0 && covid_ 
        p_covid_d_t = p_covid_d(t);
        p_covid_d_plwh_t = p_covid_d_plwh(t);
      
        % a 4x8 matrix stores # of covid death
        num_covid_d = zeros(length(p_covid_d_race), length(p_covid_d_age)); % 4x8 (race x age)
        num_covid_d_plwh = zeros(length(p_covid_d_race), length(p_covid_d_age));
        
        for i = 1:length(p_covid_d_race) % fill in the 4x8 matrix (race x age)       
            num_covid_d(i,:) = int16(sum(alive==1 & status==0) * p_covid_d_t * p_covid_d_race(i) * p_covid_d_age);          
            num_covid_d_plwh(i,:) = int16(sum(alive==1 & status~=0) * p_covid_d_plwh_t * p_covid_d_race(i) * p_covid_d_age);
        end
               
        for j = 1:length(p_covid_d_race)
            covid_ind_d = []; % stores indicies for a race
            
            for i = (1:2:length(covid_age_bucket))
                k = (i+1)/2;
                ind_d = randsample(find(alive==1 & status==0 & race==j-1 & age>=covid_age_bucket(i) &...
                            age<=covid_age_bucket(i+1)), num_covid_d(j, k)); % a vector of indecies by race by age
                
                tmp_plwh_ind = find(alive==1 & status~=0 & race==j-1 & age>=covid_age_bucket(i) &...
                            age<=covid_age_bucket(i+1));
                num = num_covid_d_plwh(j, k);        
                if length(tmp_plwh_ind) < num
                    num = length(tmp_plwh_ind);
                end   
                ind_plwh_d = randsample(tmp_plwh_ind, num); % a vector of indecies by race by age
                
                covid_ind_d = [covid_ind_d; [ind_d; ind_plwh_d]]; % a vector of indecies by race
            end
            state_matrix(covid_ind_d, StateMatCols.alive) = 0; % assign people w/ these ind alive = 0
        end         
    end
     
    %% covid setting -end-    
    
    % update policy based on year want policy to be initiated
    % this is for policies that are changing transition probabilities
    if t == policyFlag
        treatmentOn_transition_path = newPolicy_treatmentOn_transition_path;
        treatmentOff_transition_path = newPolicy_treatmentOff_transition_path;
        prepOn_transition_path = newPolicy_prepOn_transition_path;
        prepOff_transition_path = newPolicy_prepOff_transition_path;
    end


    % year tally
    yearNum = init_yr + t;
    % for now births start right away
    [state_matrix, birthTally] = birth2(state_matrix, StateMatCols, age_def, age_prop, race_def, race_prop, inflow);

    % % If we dont want birth to start right away
    % if t == 1
    %     birthTally = 0;
    % else
    %     [state_matrix, birthTally] = birth(state_matrix, StateMatCols, age_def, age_prop, race_def, race_prop, inflow);
    % end

    % Create array for case when everyone is eligible for subsetting
    all_eligible = 1:size(state_matrix, 1);

    % set a vector with enoguh values for all tally (excluding year because it is added after demog)
    yearlyTally = zeros(1, 17);
    

    % For each demographic group
    for demog_group_def = demog_table.'
        
        % demographic tally holder
        demographicTally = [];

        % Get row indices for people in the state matrix in this demographic group
        state_mat_demog_group_idx = find_demog_rows(state_matrix, StateMatCols, demog_group_def, DemogTblCols, DemogTblCols);
       
        % If there's no one in the state matrix in this demographic group, continue to the next
        if isempty(state_mat_demog_group_idx)
            continue
        end 

        %%%%%%%%%%%%%%%%%%%%%%%%%% aware %%%%%%%%%%%%%%%%%%%%%%%%%%
        %% covid testing -start-
        if ismember(t, yr_test) && test_flg %yr 2020
            % lower diagnosed probability
            aware1_transition_path = aware1_transition_path_2020;
            aware2_transition_path = aware2_transition_path_2020;
            aware3_transition_path = aware3_transition_path_2020;
        else 
            aware1_transition_path = aware1_transition_path_default;
            aware2_transition_path = aware2_transition_path_default;
            aware3_transition_path = aware3_transition_path_default;
        end
        %% covid testing -end-
        [state_matrix, awareTally1] = transition5('aware', aware1_transition_path, ...
                                    state_matrix, StateMatCols, state_mat_demog_group_idx, ...
                                    future_state_param_names, 1);
        [state_matrix, awareTally2] = transition5('aware', aware2_transition_path, ...
                                    state_matrix, StateMatCols, state_mat_demog_group_idx, ...
                                    future_state_param_names, 1);
        [state_matrix, awareTally3] = transition5('aware', aware3_transition_path, ...
                                    state_matrix, StateMatCols, state_mat_demog_group_idx, ...
                                    future_state_param_names,1);

        % this is simply to ensure that people who are aware of HIV are no longer indicated as on prep
        [state_matrix] = transition5('prep', prepClean_transition_path, ...
                                    state_matrix, StateMatCols, state_mat_demog_group_idx, ...
                                    future_state_param_names,1);

        %%%%%%%%%%%%%%%%%%%%%% tally counters %%%%%%%%%%%%%%%%%%%%%%
        
        % tally for this specific demographic group
        demographicTally = [0, 0, 0, 0, ...
                            0, 0, ...
                            awareTally1, awareTally2, awareTally3, 0, 0, 0, ...
                            0, 0, 0, 0, 0];
                            % the zero is for infection which will be added later
         
        % tally for all the demographic groups summed for this year
        yearlyTally = yearlyTally + demographicTally;
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%% ACQUIRING INFECTION %%%%%%%%%%%%%%%%%%%%%%%%%%
    %% must create a copy of the state matrix to do infections on because we do not want newly infected people to be considered in 
        % infecting new people.

    % state matrix reflecting all the changes in the year, prior to infections
    state_matrix_forYear = state_matrix;  

    % For each demographic group
    for demog_group_def = demog_table.'
        
        
        % Get row indices for people in the state matrix in this
        % demographic group
        state_mat_demog_group_idx = find_demog_rows(state_matrix, StateMatCols, demog_group_def, DemogTblCols, DemogTblCols);

       
        % If there's no one in the state matrix in this demographic group, continue to the next
        % one
        if isempty(state_mat_demog_group_idx)
            continue
        end 

        %%%%%%%%%%%%%%%%%%%%%%%%%% ACQUIRING INFECTION %%%%%%%%%%%%%%%%%%%%%%%%%%
        [state_matrix, infectionTally] = infection9(demog_group_def, DemogTblCols, mixing_matrix, mixing_table, ...
                                    MixingTblCols, state_matrix, state_matrix_forYear, StateMatCols, ...
                                    prepAdherence, prepMultiplier, numPartners, artAdherence, infectCalib, ...
                                    youngInfect, raceInfect);


        % infection is ignored when testing death rate
        % infectionTally = 0;
         % tally for this specific demographic group
        demographicTally = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, infectionTally];

        % tally for all the demographic groups summed for this year
        yearlyTally = yearlyTally + demographicTally;
    end


    % For each demographic group
    for demog_group_def = demog_table.'
        
        % demographic tally holder
        demographicTally = [];
        

        % Get row indices for people in the state matrix in this
        % demographic group
        state_mat_demog_group_idx = find_demog_rows(state_matrix, StateMatCols, demog_group_def, DemogTblCols, DemogTblCols);
   
       
        % If there's no one in the state matrix in this demographic group, continue to the next
        % one
        if isempty(state_mat_demog_group_idx)
            continue
        end 
        
        
        %%%%%%%%%%% Need to split up transions to more files %%%%%%%%%%%%%%
        %%% Examples %%%%%
        % split alive / death to death by 
            % hiv 
            % age 
            % other method (accidents)
        % split prep transition
            % prep by clinic
            % prep by hospital
            % off prep
        % treatment 
            % treatment at clinic
            % treatment at hospital
            % switch off treatment
        
        % Order is significant. How these are split and the order will depend
            % on what is found in literature and what is the naturally
            % occurence of events
        
        %%%%%%%%%%%%%%%%%%%%%%%%%% status %%%%%%%%%%%%%%%%%%%%%%%%%%
        % this modification (different from alive/death) is necessary because the change in states can
        % impact each other


        % get index of the 3 intial status that translate to the transition status

        % the identifying status is the index 1 below what it is transitioning to
            % Ex: symptomatic (2) transitions to aids (3)
        symptomaticTransIndex = find_indices(state_matrix, [1:size(state_matrix,1)], StateMatCols.status, '=', 1);
        aidsTransIndex = find_indices(state_matrix, [1:size(state_matrix,1)], StateMatCols.status, '=', 2);

        % create three new state matrix based on the three intial transitions
        symptomatic_state_matrix = state_matrix(symptomaticTransIndex,:);
        aids_state_matrix = state_matrix(aidsTransIndex,:);

        % need to create demographic group index to use because subset of the full data
        % use the find_demog_group_idx function
        symptomatic_demog_idx = find_demog_rows(symptomatic_state_matrix, StateMatCols, demog_group_def, DemogTblCols, DemogTblCols);
        aids_demog_idx = find_demog_rows(aids_state_matrix, StateMatCols, demog_group_def, DemogTblCols, DemogTblCols);

        % ONLY run the transition function IF the demog index exists
            % similar concept to the continue if the entire state matrix does not have a certain demog

        if isempty(symptomatic_demog_idx)
            statusSymptomaticTally = 0;
        else
            % update the 3 state matrix with any changes
            [symptomatic_state_matrix, statusSymptomaticTally] = transition5('status', statusSymptomatic_transition_path, ...
                                    symptomatic_state_matrix, StateMatCols, symptomatic_demog_idx, ...
                                    future_state_param_names, 1);
            % redefine the state matrix based on the index of the new matrix
            state_matrix(symptomaticTransIndex, StateMatCols.status) = symptomatic_state_matrix(:, StateMatCols.status);
        end

        if isempty(aids_demog_idx)
            statusAIDSTally = 0;
        else
            % update the 3 state matrix with any changes
            [aids_state_matrix, statusAIDSTally] = transition5('status', statusAIDS_transition_path, ...
                                    aids_state_matrix, StateMatCols, aids_demog_idx, ...
                                    future_state_param_names, 1);
            % redefine the state matrix based on the index of the new matrix
            state_matrix(aidsTransIndex, StateMatCols.status) = aids_state_matrix(:, StateMatCols.status);
        end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%% prep %%%%%%%%%%%%%%%%%%%%%%%%%%  
        % Same approach as for status

        % note that the status that transitions to prep, is not having prep (0)
        % status that transitions to no prep, is having prep (1)
        prepOnTransIndex = find_indices(state_matrix, [1:size(state_matrix,1)], StateMatCols.prep, '=', 0);
        prepOffTransIndex = find_indices(state_matrix, [1:size(state_matrix,1)], StateMatCols.prep, '=', 1);

        
        % create new state matrix based on the states
        prepOn_state_matrix = state_matrix(prepOnTransIndex,:);
        prepOff_state_matrix = state_matrix(prepOffTransIndex,:);

        % need to create demographic group index to use because subset of the full data
        % use the find_demog_group_idx function
        prepOn_demog_idx = find_demog_rows(prepOn_state_matrix, StateMatCols, demog_group_def, DemogTblCols, DemogTblCols);
        prepOff_demog_idx = find_demog_rows(prepOff_state_matrix, StateMatCols, demog_group_def, DemogTblCols, DemogTblCols);

        % ONLY run the transition function IF the demog index exists
            % similar concept to the continue if the entire state matrix does not have a certain demog
        if isempty(prepOn_demog_idx)
            prepOnTally = 0;
        else
            % update the 2 state matrix with any changes
            [prepOn_state_matrix, prepOnTally] = transition5('prep', prepOn_transition_path, ...
                                    prepOn_state_matrix, StateMatCols, prepOn_demog_idx, ...
                                    future_state_param_names, 1);
            % redefine the state matrix based on the index of the new matrix
            state_matrix(prepOnTransIndex, StateMatCols.prep) = prepOn_state_matrix(:, StateMatCols.prep);
        end

        if isempty(prepOff_demog_idx)
            prepOffTally = 0;
        else
            % update the 2 state matrix with any changes
            [prepOff_state_matrix, prepOffTally] = transition5('prep', prepOff_transition_path, ...
                                    prepOff_state_matrix, StateMatCols, prepOff_demog_idx, ...
                                    future_state_param_names, 1);
            % redefine the state matrix based on the index of the new matrix
            state_matrix(prepOffTransIndex, StateMatCols.prep) = prepOff_state_matrix(:, StateMatCols.prep);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%% treatment %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Same approach as for status

        % note that the status that transitions to treatment, is not having treatment (0)
        % status that transitions to no treatment, is having treatment (1)
        treatmentOnTransIndex = find_indices(state_matrix, [1:size(state_matrix,1)], StateMatCols.treatment, '=', 0);
        treatmentOffTransIndex = find_indices(state_matrix, [1:size(state_matrix,1)], StateMatCols.treatment, '=', 1);
        
        % create new state matrix based on the states
        treatmentOn_state_matrix = state_matrix(treatmentOnTransIndex,:);
        treatmentOff_state_matrix = state_matrix(treatmentOffTransIndex,:);

        % need to create demographic group index to use because subset of the full data
        % use the find_demog_group_idx function
        treatmentOn_demog_idx = find_demog_rows(treatmentOn_state_matrix, StateMatCols, demog_group_def, DemogTblCols, DemogTblCols);
        treatmentOff_demog_idx = find_demog_rows(treatmentOff_state_matrix, StateMatCols, demog_group_def, DemogTblCols, DemogTblCols);

        % ONLY run the transition function IF the demog index exists
            % similar concept to the continue if the entire state matrix does not have a certain demog
        if isempty(treatmentOn_demog_idx)
            treatmentOnTally = 0;
        else
            % update the 2 state matrix with any changes
            [treatmentOn_state_matrix, treatmentOnTally] = transition5('treatment', treatmentOn_transition_path, ...
                                    treatmentOn_state_matrix, StateMatCols, treatmentOn_demog_idx, ...
                                    future_state_param_names, 1);
            % redefine the state matrix based on the index of the new matrix
            state_matrix(treatmentOnTransIndex, StateMatCols.treatment) = treatmentOn_state_matrix(:, StateMatCols.treatment);
        end

        if isempty(treatmentOff_demog_idx)
            treatmentOffTally = 0;
        else
            % update the 2 state matrix with any changes
            [treatmentOff_state_matrix, treatmentOffTally] = transition5('treatment', treatmentOff_transition_path, ...
                                    treatmentOff_state_matrix, StateMatCols, treatmentOff_demog_idx, ...
                                    future_state_param_names,1);
            % redefine the state matrix based on the index of the new matrix
            state_matrix(treatmentOffTransIndex, StateMatCols.treatment) = treatmentOff_state_matrix(:, StateMatCols.treatment);
        end
        

        %%%%%%%%%%%%%%%%%%%%%% tally counters %%%%%%%%%%%%%%%%%%%%%%
        
        % tally for this specific demographic group
        demographicTally = [0, 0, 0, 0, ...
                            statusSymptomaticTally, statusAIDSTally, ...
                            0, 0, 0, prepOnTally, prepOffTally, 0, ...
                            treatmentOnTally, treatmentOffTally, 0, 0, 0];
                            % the zero is for infection which will be added later
         
        % tally for all the demographic groups summed for this year
        yearlyTally = yearlyTally + demographicTally;
    end   


    %%%%%%%%%%%%%%%%%%%%%%%%%% Add policy if flag identified %%%%%%%%%%%%%%%%%%%%%%%%%%

    % first policy is PrEP
    if policyFlag_type2 == 1

        if t >= policyFlag

            % run our intervention policy
            [prepInt_state_matrix, interventionTally] = prepCountIntervention(prepCountPolicy_transition_path, state_matrix, StateMatCols);
            
            % redefine the state matrix based on the index of the new matrix
            state_matrix(:, StateMatCols.prep) = prepInt_state_matrix(:, StateMatCols.prep);

            % add the intervention tally
            morePrepTally = [0, 0, 0, 0, ...
                            0, 0, ...
                            0, 0, 0, interventionTally, 0, interventionTally, ...
                            0, 0, 0, 0, 0];
        
            yearlyTally = yearlyTally + morePrepTally;
        end
    end


    %second policy is ART
    if policyFlag_type2 == 1

        if t >= policyFlag

            % run our intervention policy
            [artInt_state_matrix, interventionTally] = artCountIntervention(artCountPolicy_transition_path, state_matrix, StateMatCols);
            
            % redefine the state matrix based on the index of the new matrix
            state_matrix(:, StateMatCols.treatment) = artInt_state_matrix(:, StateMatCols.treatment);

            % add the intervention tally
            moreARTTally = [0, 0, 0, 0, ...
                            0, 0, ...
                            0, 0, 0, 0, 0, 0, ...
                            interventionTally, 0, interventionTally, 0, 0];
        
            yearlyTally = yearlyTally + moreARTTally;
        end
    end


    % third policy is diagnosis
    if policyFlag_type2 == 1

        if t >= policyFlag

            % run our intervention policy
            [diagInt_state_matrix, interventionTally] = diagnoseCountIntervention(diagCountPolicy_transition_path, state_matrix, StateMatCols);
            
            % redefine the state matrix based on the index of the new matrix
            state_matrix(:, StateMatCols.aware) = diagInt_state_matrix(:, StateMatCols.aware);

            % add the intervention tally
            % not that for this diag tally, these are seperate from the regularly occuring new diagnosis, we add the intervention quantity to the total
            moreDiagTally = [0, 0, 0, 0, ...
                            0, 0, ...
                            0, 0, 0, 0, 0, 0, ...
                            0, 0, 0, interventionTally, 0];
        
            yearlyTally = yearlyTally + moreDiagTally;
        end
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%% AGING %%%%%%%%%%%%%%%%%%%%%%%%%%
    % After all changes, everyone alive ages by 1year

    % Find everyone who is alive
    alive = find_indices(state_matrix, all_eligible, StateMatCols.alive, '=', 1);

    % Add a unit of time to their ages
    state_matrix(alive,StateMatCols.age) = state_matrix(alive,StateMatCols.age) + 1;
    
    % additionally, everyone with AIDS now has had aids for 1 year longer
    aidsStatus = find_indices(state_matrix, all_eligible, StateMatCols.status, '=', 3);

    % iterate if alive and aids
    alive_aids = intersect(alive, aidsStatus);

    state_matrix(alive_aids,StateMatCols.aidsYears) =  state_matrix(alive_aids,StateMatCols.aidsYears)+1;


    %%%%%%%%%%%%%%%%%%%%%%%%%% alive / death %%%%%%%%%%%%%%%%%%%%%%%%%%
    % do this at end of all other changes. want deaths to happen to population after everything,
    % For each demographic group
    for demog_group_def = demog_table.'
        
        
        % Get row indices for people in the state matrix in this
        % demographic group
        state_mat_demog_group_idx = find_demog_rows(state_matrix, StateMatCols, demog_group_def, DemogTblCols, DemogTblCols);

       
        % If there's no one in the state matrix in this demographic group, continue to the next
        % one
        if isempty(state_mat_demog_group_idx)
            continue
        end 

        %%%%%%%%%%%%%%%%%%%%%%%%%% alive / death %%%%%%%%%%%%%%%%%%%%%%%%%%
        [state_matrix, deathNaturalTally] = transition5('alive', deathNatural_transition_path, ...
                                        state_matrix, StateMatCols, state_mat_demog_group_idx, ...
                                        future_state_param_names, 1);

        % when testing aids death
        % deathNaturalTally = 0;

        [state_matrix, deathAIDSTally_aware_noTreat] = transition5_ageSplit('alive', deathAIDS_aware_noTreat_transition_path, ...
                                        state_matrix, StateMatCols, state_mat_demog_group_idx, ...
                                        future_state_param_names, deathCalib_a1, deathCalib_a2, deathCalib_a3, deathCalib_a4);
        
        % deathAIDSTally_aware_noTreat = 0;

        [state_matrix, deathAIDSTally_aware_Treat] = transition5_ageSplit('alive', deathAIDS_aware_Treat_transition_path, ...
                                        state_matrix, StateMatCols, state_mat_demog_group_idx, ...
                                        future_state_param_names, deathCalib_a1, deathCalib_a2, deathCalib_a3, deathCalib_a4);
        
        % deathAIDSTally_aware_Treat = 0;

        [state_matrix, deathAIDSTally_unaware] = transition5_ageSplit('alive', deathAIDS_unaware_transition_path, ...
                                        state_matrix, StateMatCols, state_mat_demog_group_idx, ...
                                        future_state_param_names, deathCalib_a1, deathCalib_a2, deathCalib_a3, deathCalib_a4);
        
        % deathAIDSTally_unaware = 0;

        % infection is ignored when testing death rate
        % infectionTally = 0;
         % tally for this specific demographic group
        demographicTally = [deathNaturalTally, deathAIDSTally_aware_noTreat, deathAIDSTally_aware_Treat, deathAIDSTally_unaware, ...
                             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

        % tally for all the demographic groups summed for this year
        yearlyTally = yearlyTally + demographicTally;


    end


    %%%%%%%%%%%%%%%%%%%%%%%%%% SAVE STATE MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%
    matrix_name = int2str(t);
    cd(sim_dataDir)
    save(matrix_name, 'state_matrix');
    cd(code_loc)
    
    % add births to yearly tally
    yearlyTally = [yearNum, birthTally, yearlyTally];
    
    % add row of year tallly to final tally that will be stored     
    finalTransitionTally =  [finalTransitionTally; yearlyTally];

    %%%%%%%%%%%%%%%%%%%%%%%%%% timer %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    completion = t/T*100;
    display(strcat(num2str(completion),'%'))


end

cd(sim_dataDir)
%%%%%%%%%%%%%%%%%%%%%%%%%% SAVE Tally %%%%%%%%%%%%%%%%%%%%%%%%%%
matrix_name = 'Final_Tally';
save(matrix_name, 'finalTransitionTally');
tallyHeaders = {'Year', 'Births', 'Natural Death', 'AIDS Death Diagnosed (aware) No Treatment', 'AIDS Death Diagnosed (aware) Treatment', ...
                'AIDS Death Undiagnosed (unaware)', 'Transition to Symptomatic', ...
                'Transition to AIDS', 'Transition to Aware (1)', ...
                'Transition to Aware (2)', 'Transition to Aware (3)', 'Transition On PrEP', ...
                'Transition Off PrEp', 'PrEP from Intervention', 'Transition on Treatment', 'Transition off Treatment', ...
                'Treatment from intervantion', 'Diagnosed from intervantion', 'Infection'};
csvwrite_with_headers(strcat(matrix_name,'.csv'), finalTransitionTally, tallyHeaders)

%%%%%%%%%%%%%%%%%%%%%%%%%% Counts for state matrix data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this is run rate data


% go to folder with data
% cd(dataDir)

% create list of all .mat files
% one will be Transition Tally (The csv file name for this table is Final_Tally)
% The rest are numbers reflecting the year

matlabFiles = dir('*.mat');


% load all the data files
for i = 1:length(matlabFiles)-1
    % to apply ordering, since file names start at 0, we know the names for this portion will be numeric
    % runData(i) = load(sprintf(matlabFiles(i).name));
    runData(i) = load(strcat(num2str(i-1),".mat"));
end

stateDetails = [];
for i = 1:length(matlabFiles)-1
    % get the year. Initialized with 2013 data
    % the first indexed one is year zero
    runYear = init_yr + i-1;
    % get run data
    [yearHeaders, yearData] = stateDetailsFunction(runData(i).state_matrix, StateMatCols);
    % generate desired death data and birth
    % note that in the 1st year (0 year, nobody is dead)
    if i == 1
        % overal deaths
        NumDeaths = yearData(end);
        % births
        Births = 0;
    else
        % the NumDeaths will be added to the end of state details
        % currenty cumulive number of deaths - previous cumulative number of deaths
        NumDeaths = yearData(end) - stateDetails(i-1, end-1);
        % birth factor will be added to the beginning (after run year) thus shifitng where EOY population size by 2 spts in state details
        % birth = current year death plus end of year population - end of year population from prior year
        Births = NumDeaths + yearData(1) - stateDetails(i-1, 3);
    end
    % add new row to state details
    stateDetails = [stateDetails; [runYear, Births, yearData, NumDeaths]];

    if i == length(matlabFiles)-1
        stateHeaders = ['Year', 'Births', yearHeaders, 'Number Deaths'];
    end
end

% include path to state_matrices files such that it is saved into the state matrices folder
matrix_name = 'runDetails';
save(matrix_name, 'stateDetails');
csvwrite_with_headers(strcat(matrix_name,'.csv'), stateDetails, stateHeaders)

cd(code_loc)

% toc



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% supporting function identifying the number of people in every possible compartment at every time stamp

function [outputHeaders, outputValues] = stateDetailsFunction(stateMatrix, StateMatCols)
    outputHeaders = [];
    outputValues = [];

    % susceptibel
    Susceptible = sum(stateMatrix(:,StateMatCols.alive) == 1 & stateMatrix(:,StateMatCols.status) == 0);
    % Early
    Early = sum(stateMatrix(:,StateMatCols.alive) == 1 & stateMatrix(:,StateMatCols.status) == 1);
    % Symptomatic
    Symptomatic = sum(stateMatrix(:,StateMatCols.alive) == 1 & stateMatrix(:,StateMatCols.status) == 2);
    % AIDS
    AIDS = sum(stateMatrix(:,StateMatCols.alive) == 1 & stateMatrix(:,StateMatCols.status) == 3);

    % susceptible no prep
    SuNp = sum(stateMatrix(:,StateMatCols.alive) == 1 & stateMatrix(:,StateMatCols.status) == 0 & ...
                stateMatrix(:,StateMatCols.prep) == 0);
    % susceptible prep
    SuP = sum(stateMatrix(:,StateMatCols.alive) == 1 & stateMatrix(:,StateMatCols.status) == 0 & ...
                stateMatrix(:,StateMatCols.prep) == 1);
    % Early no prep -> implies unaware
    ENp = sum( stateMatrix(:,StateMatCols.alive) == 1 & stateMatrix(:,StateMatCols.status) == 1 & ...
                stateMatrix(:,StateMatCols.prep) == 0 & stateMatrix(:,StateMatCols.aware) == 0);
    % Early prep -> implies unaware
    EP = sum(stateMatrix(:,StateMatCols.alive) == 1 & stateMatrix(:,StateMatCols.status) == 1 & ...
                stateMatrix(:,StateMatCols.prep) == 1 & stateMatrix(:,StateMatCols.aware) == 0);
    % Early no treatment -> implies aware
    ENt = sum(stateMatrix(:,StateMatCols.alive) == 1 & stateMatrix(:,StateMatCols.status) == 1 & ...
                stateMatrix(:,StateMatCols.treatment) == 0 & stateMatrix(:,StateMatCols.aware) == 1);
    % Early treatment -> implies aware
    ET = sum(stateMatrix(:,StateMatCols.alive) == 1 & stateMatrix(:,StateMatCols.status) == 1 & ...
                stateMatrix(:,StateMatCols.treatment) == 1 & stateMatrix(:,StateMatCols.aware) == 1);
    % Symptomatic no prep -> implies unaware
    SyNp = sum(stateMatrix(:,StateMatCols.alive) == 1 & stateMatrix(:,StateMatCols.status) == 2 & ...
                stateMatrix(:,StateMatCols.prep) == 0 & stateMatrix(:,StateMatCols.aware) == 0);
    % Symptomatic prep -> implies unaware
    SyP = sum(stateMatrix(:,StateMatCols.alive) == 1 & stateMatrix(:,StateMatCols.status) == 2 & ...
                stateMatrix(:,StateMatCols.prep) == 1 & stateMatrix(:,StateMatCols.aware) == 0);
    % Symptomatic no treatment -> implies aware
    SyNt = sum(stateMatrix(:,StateMatCols.alive) == 1 & stateMatrix(:,StateMatCols.status) == 2 & ...
                stateMatrix(:,StateMatCols.treatment) == 0 & stateMatrix(:,StateMatCols.aware) == 1);
    % Symptomatic treatment -> implies aware
    SyT = sum(stateMatrix(:,StateMatCols.alive) == 1 & stateMatrix(:,StateMatCols.status) == 2 & ...
                stateMatrix(:,StateMatCols.treatment) == 1 & stateMatrix(:,StateMatCols.aware) == 1);
    % AIDS no prep -> implies unaware
    ANp = sum(stateMatrix(:,StateMatCols.alive) == 1 & stateMatrix(:,StateMatCols.status) == 3 & ...
                stateMatrix(:,StateMatCols.prep) == 0 & stateMatrix(:,StateMatCols.aware) == 0);
    % AIDS Prep -> implies unaware
    AP = sum(stateMatrix(:,StateMatCols.alive) == 1 & stateMatrix(:,StateMatCols.status) == 3 & ...
                stateMatrix(:,StateMatCols.prep) == 1 & stateMatrix(:,StateMatCols.aware) == 0);
    % AIDS no treatment -> implies aware
    ANt = sum(stateMatrix(:,StateMatCols.alive) == 1 & stateMatrix(:,StateMatCols.status) == 3 & ...
                stateMatrix(:,StateMatCols.treatment) == 0 & stateMatrix(:,StateMatCols.aware) == 1);
    % AIDS treatment -> implies aware
    AT = sum(stateMatrix(:,StateMatCols.alive) == 1 & stateMatrix(:,StateMatCols.status) == 3 & ...
                stateMatrix(:,StateMatCols.treatment) == 1 & stateMatrix(:,StateMatCols.aware) == 1);
    % Cumulative Death
    CumDeath = sum(stateMatrix(:,StateMatCols.alive) == 0);

    % additional outputs of interests

    % people living with HIV
    PLWH = ENp + EP + ENt + ET + SyNp + SyNt + SyT + SyP;
    % people living with AID
    PLWA = ANp + ANt + AT + AP;
    % peoplve living with HIV/AIDS 
    PLWHA = PLWH + PLWA;
    % Population
    EOY_PopSize = sum(stateMatrix(:,StateMatCols.alive) == 1);

    % deaths must be managed outside of this function because number marked as death accumpulate
    % use the cumulative death values. The detla between years will be the desired value
    % Diagnosed deaths
    % overal deaths
    % treatment deaths

    outputHeaders = {'Population Size', 'Susceptible', 'Early', 'Symptomatic', 'AIDS', 'PLWH', 'PLWA', 'PLWHA', ...
                     'Susceptible No Prep', 'Susceptible Prep', 'Early No Prep (Undiagnosed)', ...
                     'Early Prep (Undiagnosed)', 'Early No Treatment (Diagnosed)', 'Early Treatment', 'Symptomatic No Prep (Undiagnosed)', 'Symptomatic Prep (Undiagnosed)', ...
                     'Symptomatic No Treatment (Diagnosed)', 'Symptomatic Treatment', 'AIDS No Prep (Undiagnosed)', 'AIDS Prep (Undiagnosed)', 'AIDS No Treatment (Diagnosed)', ...
                     'AIDS Treatment', 'Cumulative Death'};
    outputValues = [EOY_PopSize, Susceptible, Early, Symptomatic, AIDS, PLWH, PLWA, PLWHA, SuNp, SuP, ...
                    ENp, EP, ENt, ET, SyNp, SyP, SyNt, SyT, ANp, AP, ANt, AT, CumDeath];




end