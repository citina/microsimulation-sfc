format long
clear all
close all
clc
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to be able to run multiple iterations of the simulation in 1 go, and different simulation inputs. 
% The results will need to be saved in a unique folder hierarchy as well depending on the nubmer of iterations being run and the number of simulations being tested

% Pseudo Code
% 0) set appropriate folder
% 1) Create lists of input files to be using
% 2) Create list of the output locations (the output locations should be same length as input files list)
% 3) these should conisder running multiple iterations for a single input
% 4) run the modified simulation file 
% 	Mod file does not have any of hte analysis in it 
% 	It also dynamically pulls the input data and changes the output save destination after each iterations
% 5) run analysis iteratively over where all  the state matrices are stored. 
% 6) run cumulative analysis that averages over the multiple runs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set default values
interventionStart = 9; % Year when intervention starts 2021
iterations = (1:5);

% Flag for policy 2 starting
runSimulations = 1;
interventionFlag2 = 1;

% defining the different items being simulated
inputHeader = "Inputs";
PolicyID_p1 = "base"; %"D0" "D100"
PolicyID_p2 = "case"; %"P0_V0" "P500_V0" "P1000_V0" "P2000_V0" 
endHeader = "SF"; 
fileFormat = ".xlsx";

% defining the input file names
sim_initialPop = "init_pop_SF_V5.csv"; %%%change here

inputFiles = strings(1, length(PolicyID_p1)*length(PolicyID_p2));
tmp = reshape(1:length(inputFiles), length(PolicyID_p1), length(PolicyID_p2));
for i = 1:length(PolicyID_p1)
	for j = 1:length(PolicyID_p2)
		newFile = strcat(inputHeader, "_", PolicyID_p1(i), "_", PolicyID_p2(j), "_",endHeader,fileFormat);
		inputFiles(tmp(i, j)) = newFile; %change, was inputFiles = [inputFiles, newFile];
	end
end

% create list of output locations
workingDir = pwd;
testVersion = "SF_V3.8"; %%%change here SF_intervention
dataDirHeader = fileparts(workingDir) + "/MonteCarloResults" + "/" + testVersion; %%%change here
mkdir(dataDirHeader)

dataDirs = strings(1, length(PolicyID_p1)*length(PolicyID_p2)*length(iterations));
tmp = reshape(1:length(dataDirs), length(PolicyID_p1), length(iterations), length(PolicyID_p2));
for i = 1:length(PolicyID_p1)
	for j = 1:length(PolicyID_p2)
		for k = iterations(1):iterations(end) % one file for each iteration
			newDir = strcat(dataDirHeader, "/", PolicyID_p1(i), "/", PolicyID_p2(j), num2str(k),"/");
			dataDirs(tmp(i, k, j)) = newDir; % change, dataDirs = [dataDirs, newDir];
		end

	end
end
% dataDirs = dataDirs(dataDirs~="");

if runSimulations == 1
    
	disp("*** Running Simulation ***")
	tStart = tic; % start timer
    
	% iterate through all the input files, and then all the iterations
    for fileNum = 1:length(inputFiles)

        sim_inputFile = convertStringsToChars(inputFiles(fileNum));
        disp(strcat("*** Policy: ", num2str(fileNum), "/", num2str(length(inputFiles)), " ***"))
        disp(sim_inputFile)

        for iters = iterations(1):iterations(end)

            dataDirIndex = (fileNum-1)*length(iterations) + find(iterations == iters);
            sim_dataDir = dataDirs(dataDirIndex);
            mkdir(sim_dataDir)
            sim_dataDir = strcat(sim_dataDir, "state_matrices/"); % final output location
            mkdir(sim_dataDir)

            rng('shuffle')
            disp(strcat("*** Iteration Number: ", num2str(find(iterations == iters)),"/", num2str(length(iterations)), " ***"))
            hiv_simulation_shellMod_SF_V1; % run simulation and create state matricies

        end

    end

	
	tEnd = toc(tStart); % end timer 

	fprintf('Simulation time %d hours and %d minutes\n', floor(tEnd/3600), floor(rem(tEnd,3600)/60));

end



