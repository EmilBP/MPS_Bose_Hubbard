clear all; close all; clc;

readDirectory = '../../mnt/LinSigSeedN5fix2/';
writeDirectory = 'Plots/LinSigSeedN5fix2/';

% processStatus(readDirectory,writeDirectory);

%processBestSolutions(readDirectory,writeDirectory)
processFidelityData(readDirectory,writeDirectory);


% processCacheData(readDirectory,writeDirectory, 5.6);
