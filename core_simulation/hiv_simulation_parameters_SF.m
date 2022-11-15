
% inputWorkbook = 'Inputs_base_case_SF_V3.xlsx';

% Number of years to run the simulation
T = xlsread(inputWorkbook, 'SimDuration');

% inflow proportion
inflow = xlsread(inputWorkbook, 'InflowProportion');
% set to zero when testing  mortality rates
% inflow = 0;

% age def and proportions
age_def = xlsread(inputWorkbook, 'AgeDef')';
agePropInput = xlsread(inputWorkbook, 'AgeProp');
age_prop = agePropInput' / sum(agePropInput);

% race def and prop
race_def = xlsread(inputWorkbook, 'RaceDef')';
race_prop = xlsread(inputWorkbook, 'RaceProp')';

% Infection Calibration Constant
infectCalib = xlsread(inputWorkbook, 'InfectCalib');

% Death Calibration Constant
deathCalib = xlsread(inputWorkbook, 'DeathCalib');
deathCalib_a1 = deathCalib(1);
deathCalib_a2 = deathCalib(2);
deathCalib_a3 = deathCalib(3);
deathCalib_a4 = deathCalib(4);

% PrEP Adherence levels (not detectable, low adherence, high aherence)
% .1 multiplier means 90% efficacy. 
% 1 multiplier means 0 % efficacy
prepAdherence = xlsread(inputWorkbook, 'PrepAdherence')';

% PrEP infection multipler 
prepMultiplier = xlsread(inputWorkbook, 'PrepMult')';

% ART Adherence
artAdherence = xlsread(inputWorkbook, 'ArtAdherence');

% Multiplier for increased force of infection for young
youngInfect = xlsread(inputWorkbook, 'YoungInfect');

% Multiplier for increased force of infection based on race
raceInfect = xlsread(inputWorkbook, 'RaceInfect');

[~, filePaths] = xlsread(inputWorkbook, 'FilePaths', 'B1:B37');

% Initial population state matrix input file
file_common_part = fileparts(code_loc) + "/";
init_pop_file = file_common_part + filePaths{1};

% File path for table that defines all demographic variable categories
demog_var_def_file = file_common_part + filePaths{2};
% for testing death rate
%demog_var_def_file = 'input/demog_var_def_testingDeathRate.csv';

% File path for table that defines all the demographic variables for mixing
% matrix
mixing_mat_def_file = file_common_part + filePaths{3};

% File path for mixing matrix
mixing_mat_file = file_common_part + filePaths{4};

% File path for number of partners based on mixing stratification;
num_partners_file = file_common_part + filePaths{5};

% Folder path for where the state matrices are saved
state_matrices_path = file_common_part + filePaths{6};

% Paths defined for aware transitions
aware1_transition_path = file_common_part + filePaths{7};
aware2_transition_path = file_common_part + filePaths{8};
aware3_transition_path = file_common_part + filePaths{9};
% lower diagnosed probability under covid 
% if isempty(filePaths{35})==0
try 
    aware1_transition_path_2020 = file_common_part + filePaths{35};
    aware2_transition_path_2020 = file_common_part + filePaths{36};
    aware3_transition_path_2020 = file_common_part + filePaths{37};
catch
end

% Paths defined for death transition
deathAIDS_aware_noTreat_transition_path = file_common_part + filePaths{10};
deathAIDS_aware_Treat_transition_path = file_common_part + filePaths{11};
deathAIDS_unaware_transition_path = file_common_part + filePaths{12};
deathNatural_transition_path = file_common_part + filePaths{13};

% Paths defined for prep transition
prepOn_transition_path = file_common_part + filePaths{14};
prepOn_transition_path_2013 = file_common_part + filePaths{15};
prepOn_transition_path_2014 = file_common_part + filePaths{16};
prepOn_transition_path_2015 = file_common_part + filePaths{17};
prepOn_transition_path_2016 = file_common_part + filePaths{18};
prepOn_transition_path_2017 = file_common_part + filePaths{19};
prepOn_transition_path_2018 = file_common_part + filePaths{20};
prepOff_transition_path = file_common_part + filePaths{21};
prepClean_transition_path = file_common_part + filePaths{22};

% Paths defined for status
statusAIDS_transition_path = file_common_part + filePaths{23};
statusSymptomatic_transition_path = file_common_part + filePaths{24};

% Paths defined for treatment
treatmentOn_transition_path = file_common_part + filePaths{25};
treatmentOff_transition_path = file_common_part + filePaths{26};
aware3_treatment_diagnosis_transition = file_common_part + filePaths{27};

% Paths for policy
newPolicy_treatmentOn_transition_path = file_common_part + filePaths{28};
newPolicy_treatmentOff_transition_path = file_common_part + filePaths{29};
prepCountPolicy_transition_path = file_common_part + filePaths{30};
newPolicy_prepOn_transition_path = file_common_part + filePaths{31};
newPolicy_prepOff_transition_path = file_common_part + filePaths{32};
artCountPolicy_transition_path = file_common_part + filePaths{33};
diagCountPolicy_transition_path = file_common_part + filePaths{34};

% Define future state parameter names
future_state_param_names = {'state', 'prob'};










