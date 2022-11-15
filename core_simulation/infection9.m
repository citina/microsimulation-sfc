% This function will take in a particular demographic group, a state
% matrix, and a mixing matrix and will transition people in this
% demographic group to infection based on a uniform random number draw

% infection requries 2 state matrix. One is the state matrix being updated with state progression. (the output)
% the 2nd state matrix is the one that represents the status at the beginning of the time window (year). 
    % this is necessary because if the state matrix reflecting progression is being used ot calculated probability of infection,
    % the chance of infection is higher for each sequential demographic group, because those who were only susceptible before have become infected.

function [updated_state_matrix, infectionTally] = infection8(reference_demog_group_def, DemogTblCols, mixing_matrix, mixing_table, ...
                                                             MixingTblCols, state_matrix_toUpdate, state_matrix_forYear, StateMatCols, prepAdherence, ...
                                                             prepMultiplier, numPartners, artAdherence, calibrationConstant, ...
                                                             YoungForceOfInfectionScale, RaceForceOfInfectionScale)
    
    % make a copy of the state matrix such that a tally can be made of the number of changes that occur
    state_matrix_copy = state_matrix_toUpdate;

    % The reference dograpchi group is formatted as [age min, age max, race].
    % we can thus check the 2nd number and if it is under 25, we know to divide by our risk due to condom use
    % and apply it to our calibration constant... 
    % this works because our demographic group cutoff is at 25

    % if young, divide by the scale (becuase scale less than one)
    if reference_demog_group_def(2) < 25
        ForceOfInfection = calibrationConstant / YoungForceOfInfectionScale;
    else
        ForceOfInfection = calibrationConstant;
    end
    % if not young, leave calibration constant as is


    % the force of infection must then be scaled by the race characteristic.
    % the race characteristic is the third one of demographic reference (must add one cuz demog group starts at 0)
    ForceOfInfection = ForceOfInfection * RaceForceOfInfectionScale(reference_demog_group_def(3) + 1);



    % Returns the row index based on the mixing_mat_def
    % demographic group table
    demog_row = find_demog_rows(mixing_table, MixingTblCols, reference_demog_group_def, DemogTblCols, MixingTblCols);

    % Based on the demographic row, extract the mixing probabilities from the
    % corresponding row number in the mixing matrix
    mixingProbabilites = mixing_matrix(demog_row,:);

        
    % Placeholder for accounting for all the infection probabilities over all
    % demographic group our reference group is mixing with
    total_no_infec_prob = 1;
    
    % For each of all possible demographic groups
    for demog_group_idx = 1:size(mixingProbabilites, 2)

        % Get the demographic group definition we are mixing with
        mix_demog_group_def = mixing_table(demog_group_idx,:);
        % pull mixing probability for the particular instnace (who the susceptible person is mixing with)
        mixingProb_group = mixingProbabilites(demog_group_idx);
        % pull number of partners associated with demographic row (the demographic of suceptible person)
        numPart = numPartners(demog_row);

        % determine probability of no infection when mixing with this particular race
        % account for the number of partners in the particular racial group
        % calibration constnat is applied to the I/N calculation (calf_infec_prob)
        no_infec_prob = (1 - ForceOfInfection*calc_infec_prob4(state_matrix_forYear, StateMatCols, mix_demog_group_def, ...
                                             MixingTblCols, MixingTblCols, artAdherence)) ^ (numPart * mixingProb_group);
        
        % Multiply this to the total no infection probability over all stratifications may mix with
        % groups we are mixing with
        total_no_infec_prob = total_no_infec_prob * no_infec_prob;
    
    end  

  
    % Find the people in the demographic group
    reference_demog_mat_indices = find_demog_rows(state_matrix_toUpdate, StateMatCols, reference_demog_group_def, DemogTblCols, MixingTblCols);
    
    % Find eligible people to infect in our reference demographic group
    % Susceptible individuals have: Status susceptible, alive
    susceptible = find_indices(state_matrix_toUpdate, reference_demog_mat_indices, StateMatCols.alive, '=', 1);
    susceptible = find_indices(state_matrix_toUpdate, susceptible, StateMatCols.status, '=', 0);


    % PrEP is a multiplier considered towards infection
    % as a multiplier, it is assumed the prep has equal efficacy across all demographics
    % differentiate susceptible with prep vs susceptible without prep
    susceptible_prep = find_indices(state_matrix_toUpdate, susceptible, StateMatCols.prep, '=', 1);
    susceptible_noprep = find_indices(state_matrix_toUpdate, susceptible, StateMatCols.prep, '=', 0);




    % infected without prep

    % Do a random number draw for all of the eligible rows not on prep
    unif_prob = rand([1, length(susceptible_noprep)]);

    %prob of infection for individual in this demographic group and no prep aspect
    probInfection = 1 - total_no_infec_prob;



    % Compare to random numbers picked from a uniform
    % distribution and determine who is infected
    to_infect = probInfection > unif_prob;
                        
    % Get indices to infect based on people who are eligible
    state_mat_infec_rows = susceptible_noprep(to_infect);
    
    % Infect chosen people
    state_matrix_toUpdate(state_mat_infec_rows, StateMatCols.status) = 1;


    % just to see if prob infection changing
    % probInfection


    % infected with prep

    % need to create a vector that shows the prep multiplier levels determined by user adherence 
    % create a for loop to iterate through the levels of prep adherence
    adherencePrepMultArray = [];
    for i = 1:length(prepAdherence)
        adherence = ones(1, round(length(susceptible_prep)*prepAdherence(i)))*prepMultiplier(i);
        adherencePrepMultArray = [adherencePrepMultArray; adherence'];
    end

    % check to make sure length of adherence multArray is approporiate
    adherencePrepMultArray = checkLength(adherencePrepMultArray, length(susceptible_prep), prepMultiplier);

    % randomize order of prep multiplier
    adherencePrepMultArray = adherencePrepMultArray(randperm(numel(adherencePrepMultArray)));

    % Do a random number draw for all of the eligible rows on prep
    unif_prob = rand([1, length(susceptible_prep)]);

    % we assume all individuals have same average number of partners
    % prob to infect when have multiple partners is binomial with p = probability infection by 1 person
    % we want 1 - p(infection by nobody) to show probability that infection by at least 1 person
        % this could need to be swithced to a distribution such that each person has dif number of partners

    %prob of infection for individual in this demographic group and no prep aspect
    probInfection = (1 - total_no_infec_prob)*adherencePrepMultArray;


    % Compare to random numbers picked from a uniform
    % distribution and determine who is infected
    to_infect = probInfection > unif_prob';
                        
    % Get indices to infect based on people who are eligible
    state_mat_infec_rows = susceptible_prep(to_infect);
    
    % Infect chosen people
    state_matrix_toUpdate(state_mat_infec_rows, StateMatCols.status) = 1;

    
    % Return the updated state matrix
    updated_state_matrix = state_matrix_toUpdate;

    % return infection infection infectionTally

    infectionTally = sum(state_matrix_copy(:, StateMatCols.status) ~= state_matrix_toUpdate(:, StateMatCols.status));

    % display(reference_demog_group_def)    
    % infectionTally


    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% supporting function


function outputVect = checkLength(inputVector, desiredLength, supValArray)
    outputVect = inputVector;
    if length(inputVector) ~= desiredLength
        % if too many, remove last one
        if length(inputVector) > desiredLength
            outputVect = inputVector(1:desiredLength,1);
        % if too few, add one
        else
            outputVect = [inputVector; (randsample(supValArray, desiredLength - length(inputVector), true))'];
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
