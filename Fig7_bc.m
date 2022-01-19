% Fig7_bc
% Code to generate MUA and PSTH example data from optogenetic stimulation of V6

% load MUA data from V1 and V6
addpath data
load('Fig7_bc_data.mat')

% plotting parameters
preT = 0.6;     % pre-stimulus time
postT = 1.25;   % post-stimulus time
Fs = 1000;      % sampling rate
xlimPl = [-0.1 0.25+0.1]; % xaxis for plotting
msSmooth = 2;   % smoothing window size
plotName  = 'Fig7bc';

% reject noisy trials with large artifacts based on std
cutoffSTD = 10;
cleanSTDtrials = trlSTDSelec(aMU1eS, aMU6eS, cutoffSTD);

% plot MUA PSTH
[stats] = plotMUAV1V6meanOpto_ksTest(plotName, tr, aMU1eS, aMU6eS, preT, postT, xlimPl, Fs, msSmooth);
