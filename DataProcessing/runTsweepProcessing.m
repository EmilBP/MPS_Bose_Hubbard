clear all; close all; clc;

readDirectory = '../../../mnt/LinSigSeedN5fix/';
writeDirectory = '../../DataProcessing/Plots/LinSigSeedN5fix/';

% processStatus(readDirectory,writeDirectory);

%processBestSolutions(readDirectory,writeDirectory)
processFidelityData(readDirectory,writeDirectory);


% processCacheData(readDirectory,writeDirectory, 5.6);
