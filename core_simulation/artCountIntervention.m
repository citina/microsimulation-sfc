


function [updated_state_matrix, interventionTally] = artCountIntervention(interventionPath, state_matrix, StateMatCols)


    % make a copy of the state matrix such that a tally can be made of the number of changes that occur
    state_matrix_copy = state_matrix;


    % Read in intervention plan
   [intervention_mat, interventionCols] = read_table(interventionPath);


   % get index for the intervetnion matrix (not same as state matrix)
   % this is hardcoded but we could use fieldnames (see transition.m)
   ageMin_index = 1;
   ageMax_index = 2;
   % we will ignore race and status for this intervention
%    race_index = 3;
%    status_index = 4;
   aware_index = 3;
   treat_index = 4;
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
       % we ignore race and status
    %    race = subgroup(race_index);
    %    status  = subgroup(status_index);
       aware = subgroup(aware_index);
       treatment = subgroup(treat_index);
       alive = subgroup(alive_index);

       % number of ART
       numART = subgroup(quant_index);

       % get index of people intersted in
       aliveVec = state_matrix(:,StateMatCols.alive);
       ageVec = state_matrix(:,StateMatCols.age);
         %  raceVec = state_matrix(:,StateMatCols.race);
       % dont need status vect. if aware, cannot be undiagnosed
       statusVec = state_matrix(:,StateMatCols.status);
       awareVec = state_matrix(:,StateMatCols.aware);
       treatVec = state_matrix(:,StateMatCols.treatment);

       % get the index of all individuals who fit criteria
      % hard code that status must not be 0
       boolVec = (ageVec >= ageMin & ageVec <= ageMax & awareVec == aware & treatVec == treatment & aliveVec == alive);

       % get indices of interest
       indices = find(boolVec);

       length(indices);

       % if the number of people who can start ART is smaller than desired amount to give
           % set the num on prep to the min
       numART = min([numART, length(indices)]);

       % randomly sample the number of individuals without replacement
       giveART = datasample(indices, numART, 'Replace', false);

       startART = ones(numART, 1);

       % change the indices of the people identified
       state_matrix(giveART, StateMatCols.treatment) = startART;

   end


   % double check that the proper number of changes were made

   A = state_matrix_copy(:, StateMatCols.treatment) ~= state_matrix(:, StateMatCols.treatment);

   interventionTally = sum(A);

   % display("This is the number of ART we added via intervention")
   % interventionTally

   updated_state_matrix = state_matrix;

end


