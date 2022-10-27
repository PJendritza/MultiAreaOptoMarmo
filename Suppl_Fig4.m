% Code to produce Supplementary Fig. 4
% Yellow vs. blue laser control for Monkey D
% Session: 20573, sine wave stimulation

%% load data
load Suppl_Fig4_data.mat

%% plotting parameters
figName = 'Monkey D opto control';
preT = 0.6;
postT = 1.25;
Fs = 1000;
xlimPl = [-0.1 0.25+0.1];
msSmooth = 2;

[stats] = plotMUAV1V6meanOpto_ksTest_MonkeyD_control( figName, tr, aMU1eS, aMU6eS, preT, postT, xlimPl, Fs, msSmooth);

%% plot blue laser signal, scaled to mW amplitude
peak_mW = 25; % laser amplitude was 25mW in all trials
scaledLaserSignal505 = rescale(ledSignal,0,peak_mW);

timeXa = linspace(-preT,postT,size(scaledLaserSignal505,2)); % time for plotting

figure
plot(timeXa, scaledLaserSignal505, 'color', [0.0 0.62 0.72])
xlim(xlimPl)
ylim([-4 30])
ylabel('Power (mW)')
xlabel('Time (s)')
title('505 nm laser signal')