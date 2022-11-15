clear
clear all
close all
clc
format long g

addpath([fileparts(pwd), '/global_fun']);

% create vectors for state, race, and age that contain the actual number of individuals infected
% This is based on 2012, the first year of infection
% NHAS 2013
diag_state_infected = [144, 151, 50]; % >500, 200-500, <=200
diag_race_infected = [37, 87, 164, 57]; %B,H,W,O
diag_age_infected = [119, 189, 36, 1];

assert(isequal(sum(diag_state_infected), sum(diag_race_infected)));
assert(isequal(sum(diag_state_infected), sum(diag_age_infected)));
assert(isequal(sum(diag_age_infected), sum(diag_race_infected)));

% Infected and undiagnosed = Population size * prevalence * proportion undiagnosed
N = 58204; %changed
PI = 0.23; %changed
PUDI = 0.07; 
totalPop = N*PI*PUDI; %total infected and undiagnosed

% by state, by race
pop_state_infected = [.413, 0.503, 0.084]; %same as SD
pop_race_infected = [0.09, 0.19, 0.63, 0.09]; %B,H,W,O %diag plwh NHAS
pop_age_infected = [0.24, 0.41, 0.26, 0.09]; %general population, diag age is overtime


% create a diagnoal matrix that, when multiplied to our objective, results in a vector of people diangosed in each bucket

daigVector = [];

for i = 1:length(pop_state_infected)

	for j = 1:length(pop_race_infected)

		for k = 1:length(pop_age_infected)

			% our initial population assumes indepdence of these 3 criteria in our initial population for those who are undiagnosed
			daigVector = [daigVector, totalPop * pop_state_infected(i)*pop_race_infected(j)*pop_age_infected(k)];

		end
	end
end

popMatrix = diag(daigVector);



% create vectors to extract values for state, race, and age
% use sparse matrix

% states 
% the frist 16 entries are state 1, next 16 are state 2, and next 16 are state 3
i = [1:1:48];
j = [ones(1,16), ones(1,16)*2, ones(1,16)*3];
k = ones(1,48);
extract_states = sparse(i, j, k);


% race
% 1-4, 13-16, 25-28 are race 1
i = [1:4, 17:20, 33:36, 5:8, 21:24, 37:40, 9:12, 25:28, 41:44, 13:16, 29:32, 45:48];
j = [ones(1,12), ones(1,12)*2, ones(1,12)*3, ones(1,12)*4];
k = ones(1,48);
extract_race = sparse(i, j, k);


% age
% every 4th entry is the same age bucket
i = [1:4:48, 2:4:48, 3:4:48, 4:4:48];
j = [ones(1,12), ones(1,12)*2, ones(1,12)*3, ones(1,12)*4];
k = ones(1,48);
extract_age = sparse(i, j, k);


cvx_begin

	variable diagRates(1,48)
	variable state_errors(1,3)
	variable race_errors(1,4)
	variable age_errors(1,4)

	% Objective
	% minimize(sum_square(state_errors) + sum_square(race_errors) + sum_square(age_errors))

	% Objective 2
	% minimize(3/3*sum_square(state_errors) + 2/3*sum_square(race_errors) + 1/4*sum_square(age_errors))
	% state and race have 3 levels, but age has 4. Each is term is divided by the nember of levels. 
	% we want stage to be more important than race to be more importnat than age

	% objective 3
    % minimize(1/3*sum_square(state_errors) + 2/3*sum_square(race_errors) + 2/4*sum_square(age_errors))
    
    
    % objective 5 (weight such that each term has same weight [4 races, 3 stages, 4 ages] )
	minimize(1/4*sum_square(state_errors) + 1/3*sum_square(race_errors) + 1/4*sum_square(age_errors))

	subject to

		diagRates >= 0;
		diagRates <= 1;

		diagnosedCount = diagRates * popMatrix;

		diag_state_infected == diagnosedCount * extract_states + state_errors	
		diag_race_infected == diagnosedCount * extract_race + race_errors	
		diag_age_infected == diagnosedCount * extract_age + age_errors	
		

cvx_end

diagRates'

% save('diagnosedProportions', 'diagnosedRates')
cd(fileparts(pwd)+"/input/transitions")
save('diagnosedProportions_SF_V2', 'diagRates')

