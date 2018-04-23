clear all; close all; clc;

readDirectory = '../../../mnt/LinSigSeedN5fix3/';
writeDirectory = '../../DataProcessing/Plots/LinSigSeedN5fix3/';

% processStatus(readDirectory,writeDirectory);


% processBestSolutions(readDirectory,writeDirectory)

processFidelityData(readDirectory,writeDirectory);


% for T = 1:0.5:4.5
%     processCacheData(readDirectory,writeDirectory, T);
% end