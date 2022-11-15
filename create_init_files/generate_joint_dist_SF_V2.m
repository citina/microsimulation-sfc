clear
clear all
close all
clc
format long g

addpath([fileparts(pwd), '/global_fun']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this function is create the joint distribution we should be using for the initial 
% population in our simulation

% Literature typically presents populations or demographics stratified by single characteristics, and
% not necessrily joint.

% We have population data for LA County stratified by single characteristics over the span of 6 years.
% We would like to determine the joint distribution for the initial population we will use in our model.

% we will assume linear transfermations between years, to simulate disease progression over the years with available data
% We formulate the following optimization problem to identify the joint distribution of our initial population

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%% Decision Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% race_age_vls 
% row vector that contains the number of people in each stratification for the initial year. 
% This will be the inital population we are trying to determine
% Stratification levels include the following
	% Age: 15-19, 20-29, 30-39, 40-49, 50-59, 60-69, 70-79, 80-99
	% Race: Black, Hispanic, White, Other
	% Treatment / Viral Suppression: Yes, No
% total of 32 unknowns for initial population

% slack variables will also be decision variables as part of our objective
% each year will have a slack variable vector for both age (4 levels) and race (4 levels) vls constraitns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLWH
% total initial population: 13126 from NHAS

% ageBucket
% number of people in each age bucket of initial year: 

% raceBucket
% number of people in each race bucket of initial year:


% The following parameters exist for each year of interest

% vls_age
% number of people virally suppressed in each age bucket


% vls_race
% number of people virally suppressed in each race bucket (re-normalized for used races)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%% Constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% race_age_vls vector must be non-negative


% matching race and age numbers for year 1
% the intiial year is reflected via the design variables
	% all other years are effectively transformations of the intial year
% sum of race_age_vls over each age bucket must match the age proportions
% sum of race_age_vls over each race bucket must match the race proportions


% matching race and age viral suppresion proportions
% these can be strict equality with or without slack variabels
% sum of race_age_vls over each age and vls bucket must match the age and vls proportion
% sum of race_age_vls over each race and vls bucket must match the race and vls proporiton
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%% objective %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% multiple objectives to test
% minimize error at over all years
% weighted error sthat should be minimized
% forced matching on certain levels and minimizing errors on others

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%% transformations for age, race, and vls for constraints %%%%%%%%%%%%%%%%%%%
% A1-A8, R1, T1 
% A1-A8, R2, T1
% A1-A8, R3, T1 
% A1-A8, R4, T1
% A1-A8, R1, T2 
% A1-A8, R2, T2
% A1-A8, R3, T2 
% A1-A8, R4, T2

% we will use the sparse function to define all our square matrices

% each inner vector refers to one of the age buckets
i = [[1:4 33:36], [5:8 37:40], [9:12 41:44], [13:16 45:48] [17:20 49:52] [21:24 53:56] [25:28 57:60] [29:32 61:64]]; 
j = repelem([1 2 3 4 5 6 7 8], 8);
k = repelem(1, 64);
extract_age = sparse(i, j, k); 

i = [1:4:64, 2:4:64, 3:4:64, 4:4:64]; 
j = repelem([1 2 3 4], 16);
k = repelem(1, 64);
extract_race = sparse(i, j, k); 

i = [1:4, 5:8, 9:12, 13:16, 17:20, 21:24, 25:28, 29:32]; 
j = repelem([1 2 3 4 5 6 7 8], 4);
k = repelem(1, 32);
extract_vls_age = sparse(i, j, k, 64, 8); 

i = [1:4:32, 2:4:32, 3:4:32, 4:4:32]; 
j = repelem([1 2 3 4], 8);
k = repelem(1, 32);
extract_vls_race = sparse(i, j, k, 64, 4); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% RHS of contraints (parameters we compare against) %%%%%%%%%%%%%%%

% data from NHAS SF compressed 2012
PLWH_2012 = round(13126*0.95); %changed

% number of PLWH by age and by race
% change: use more detailed numbers (NHAS 2012 data)
age_2012 = round([7, 545, 1701, 4418, 4292, 1829, 303, 31]*0.95); 
race_2012 = round([0.091, 0.186, 0.633, 0.090] * PLWH_2012); %B,H,W,O

% number of people virally suppressed by age and by race
vls_age_2012 = round([0.5, 0.4977, 0.5910, 0.6370, 0.6923, 0.7508, 0.7483, 0.7] .* age_2012); %?
vls_race_2012 = round([.5864, .6577, .6768, 0.6661] .* race_2012); %B,H,W,O %?


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% normalizing age and weight factors
ageWeight_vls = 1/8;
raceWeight_vls = 1/4;

% weight for overall age and race
ageWeight = 1/8;
raceWeight = 1/4;


cvx_begin 

	% cvx_solver SeDuMi
	% cvx_solver Gurobi

	% decision variables
	variable P_2012(1, 64) 


	% need to include slack variables for race and age in 2014 and birth
	% because multiplying proportions by total PLWH, it is possible that with rounding, numbers are off 1 by 1
	% between age and race.
	% this holds true for the 2012 initial population as well
	
	variable errors_age(1, 8)
	variable errors_race(1, 4)


	% errors
	variable errors_age_vls_2012(1, 8) 
	variable errors_race_vls_2012(1, 4) 



	% objective


	% weights assinged in this manner because age has 3 stratifications while race has 4
	% includes vls and age and race over time as criteria
	minimize(ageWeight_vls*sum_square(errors_age_vls_2012) + ...
			 raceWeight_vls*sum_square(errors_race_vls_2012) + ...
			 ageWeight*sum_square(errors_age) + ...
			 raceWeight*sum_square(errors_race))

    

	% constraints
	subject to

		
		% initial year

		% non-negative
		P_2012 >= 0;


		% slack variables must be less than or equal to 1 and greater than or equal to 0
		% this only applies to our initial population, which is fixed, and our births.
		% Does NOT apply to populations that occur in later years

		% P_2014_age_slack <= ones(1,4) 
		% P_2014_age_slack >= -1 *  ones(1,4)
		% P_2014_race_slack <= ones(1,3) 
		% P_2014_race_slack >= -1 *  ones(1,3)  		

		% % initial age matches
		P_2012 * extract_age + errors_age == age_2012;
		P_2012 * extract_race + errors_race == race_2012;
 		
 		% 2012, initial year, race and age vls numbers
 		P_2012 * extract_vls_age + errors_age_vls_2012 == vls_age_2012;
 		P_2012 * extract_vls_race + errors_race_vls_2012 == vls_race_2012;



cvx_end

	

	
diagnosedPopulation = round(P_2012');

save('DiagnosedPop_2012_crossSectional_SF_V2', 'diagnosedPopulation')

csvwrite("tmp2.csv", diagnosedPopulation)

