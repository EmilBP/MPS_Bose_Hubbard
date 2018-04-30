clear all; close all; clc;

readDirectory = '../../../mnt/LinSigSeedN20Tsweep/';
writeDirectory = '../../DataProcessing/Plots/LinSigSeedN20Tsweep/';

% processStatus(readDirectory,writeDirectory);


processBestSolutions(readDirectory,writeDirectory)

processFidelityData(readDirectory,writeDirectory);


% for T = 1:0.5:4.5
%     processCacheData(readDirectory,writeDirectory, T);
% end