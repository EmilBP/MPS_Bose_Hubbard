clear all; close all; clc;

readDirectory = '../../../mnt/LinSigSeedN5fix3/';
writeDirectory = '../../DataProcessing/Plots/LinSigSeedN5fix3/';

processStatus(readDirectory,writeDirectory);

%processBestSolutions(readDirectory,writeDirectory)
processFidelityData(readDirectory,writeDirectory);


% processCacheData(readDirectory,writeDirectory, 5.6);
