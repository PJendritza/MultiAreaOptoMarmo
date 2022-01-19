% Figure Fig5_cdfg
% Code to generate MUA and PSTH example data from area V1 and V6

% load MUA data from V1 and V6
addpath data
load('Fig5_cdfg_data.mat')

% plotting parameters
preT = 0.6;     % pre-stimulus time
postT = 1.25;   % post-stimulus time
Fs = 1000;      % sampling rate
xlimPl = [-0.25 0.65+0.3]; % xaxis for plotting
msSmooth = 8;   % smoothing window size
plotName  = 'Fig5cdfg';

% reject noisy trials with large artifacts based on std
cutoffSTD = 10; % standard deviation cutoff for noisy trials
cleanSTDtrials = trlSTDSelec(aMU1eS, aMU6eS, cutoffSTD); % 

% plot MUA PSTH
[stats] = plotMUAV1V6mean_ksTest( plotName, aMU1eS(cleanSTDtrials,:,:), aMU6eS(cleanSTDtrials,:,:), preT, postT, xlimPl, Fs, msSmooth);





