format long
close all
clear all
clc
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A wrapper for all covid scenarios
% v1: all HIV care reduction on until 2022, partners reduction for year 2020
% v2: all HIV care reduction on until 2025, partners reduction for year 2020
% v3: all HIV care reduction on, partners reduction off for year 2020
% v4: all HIV care reduction off, partners reduction on for year 2020
% v7 - resources for existing patients: HIV testing and prep initiation
%      reduction on until 2025, VLS reduction and prep discontinuation low
%      adherence on until 2022, partners reduction for year 2020
% v8 - resources for new patients: HIV testing and prep initiation
%      reduction on until 2022, VLS reduction and prep discontinuation low
%      adherence on until 2025, partners reduction for year 2020
% 
% 2020:8, 2022:10, 2025:13, 2030:18, 2035:23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set default values
iterations = (1:1);
covid_ = 1;

% create list of output locations
input_path = "MonteCarloResults/" + "SF_covid_test";
dataDirHeader = pwd + "/"+ input_path;
mkdir(dataDirHeader)

% v0 - base case
[yr_test, yr_prep_init, yr_vls, yr_prep_dis, yr_par, yr_covid] = deal(8:8, 8:8, 8:8, 8:8, 8, 8:8);
[test_flg, prep_init_flg, vls_flg, prep_dis_flg, par_flg] = deal(0, 0, 0, 0, 0);

% v1
% [yr_test, yr_prep_init, yr_vls, yr_prep_dis, yr_par, yr_covid] = deal(8:10, 8:10, 8:10, 8:10, 8, 8:10);
% [test_flg, prep_init_flg, vls_flg, prep_dis_flg, par_flg] = deal(1, 1, 1, 1, 1);

% v2
% [yr_test, yr_prep_init, yr_vls, yr_prep_dis, yr_par, yr_covid] = deal(8:13, 8:13, 8:13, 8:13, 8, 8:13);
% [test_flg, prep_init_flg, vls_flg, prep_dis_flg, par_flg] = deal(1, 1, 1, 1, 1);

% v3
% [yr_test, yr_prep_init, yr_vls, yr_prep_dis, yr_par, yr_covid] = deal(8:8, 8:8, 8:8, 8:8, 8, 8:8);
% [test_flg, prep_init_flg, vls_flg, prep_dis_flg, par_flg] = deal(1, 1, 1, 1, 0);

% v4
% [yr_test, yr_prep_init, yr_vls, yr_prep_dis, yr_par, yr_covid] = deal(8:8, 8:8, 8:8, 8:8, 8, 8:8);
% [test_flg, prep_init_flg, vls_flg, prep_dis_flg, par_flg] = deal(0, 0, 0, 0, 1);

% v7 - resources for existing patients
% [yr_test, yr_prep_init, yr_vls, yr_prep_dis, yr_par, yr_covid] = deal(8:13, 8:13, 8:10, 8:10, 8, 8:13);
% [test_flg, prep_init_flg, vls_flg, prep_dis_flg, par_flg] = deal(1, 1, 1, 1, 1);

% v8 - resources for new patients
% [yr_test, yr_prep_init, yr_vls, yr_prep_dis, yr_par, yr_covid] = deal(8:10, 8:10, 8:13, 8:13, 8, 8:13);
% [test_flg, prep_init_flg, vls_flg, prep_dis_flg, par_flg] = deal(1, 1, 1, 1, 1);

% memo
[covid_til, test_til, prep_init, vls_til, prep_dis, par_til] = ...
    deal(2012+yr_covid(end), 2012+yr_test(end), 2012+yr_prep_init(end), ...
    2012+yr_vls(end), 2012+yr_prep_dis(end), 2012+yr_par(end));
tmp_T = table(covid_til, test_til, prep_init, vls_til, prep_dis, par_til, ...
    test_flg, prep_init_flg, vls_flg, prep_dis_flg, par_flg);
writetable(tmp_T, dataDirHeader+'/memo.txt', 'Delimiter', '|')

SImulationShellScript_SF_covid;

