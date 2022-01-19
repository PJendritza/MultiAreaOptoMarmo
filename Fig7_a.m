% Fig7_a
% Code to generate spiking example traces during optogentic stimulation in V6


% session parametrs 
cfg.sessionID  = 475;
figName  = 'Fig7a';

% start end end time when data was recorded
sessionStartT = 200; 
sessionEndT = 860; 

% load data
addpath data
load('Fig7_a_data.mat')

% define channel mapping
mapnames = {'V1_1A'; 'V1_1B'; 'V6_1A'; 'V6_1B'; 'V6_2A'; 'V6_2B'};

% select channels to plot
cfg.selShankNr    = 5; %  select which shank to load here (1 to 6)
chp_V6 = [16:26]; % plot these channels  from the shank

mappname_V6 = mapnames{cfg.selShankNr}; % select which shank to load here
chns_V6 = getMappingV1V6( mappname_V6 ,'raw',cfg.sessionID)';
chnsNormal_V6 = getMappingV1V6( mappname_V6 ,'normal',cfg.sessionID)';

%% define time vector
FsRaw = 25000;

% set start and end times for plotting
startRawT = 508.75;
endRawT = startRawT+1;

startRawI = (startRawT-sessionStartT)*FsRaw;
endRawI =   (endRawT-sessionStartT)*FsRaw-1;

timePlotRaw = linspace(startRawT,endRawT, endRawI-startRawI+1);

FsLEDV  = 25000;  
TimeContLEDV = linspace(1/FsLEDV,1/FsLEDV*numel(LEDV(:,1)),numel(LEDV(:,1)));

startLEDVI = (startRawT)*FsLEDV;
endLEDVI =   (endRawT)*FsLEDV-1;

%% filter data
HighPass = 300;
LowPass = 6000;

[Z,P,K] = butter(4, [HighPass LowPass]/(FsRaw/2)); % 4 or 2?
[sosF, gain] = zp2sos(Z,P,K);
disp('Filtering with 4th order Butterworth filter (300-6000Hz)')

dataV6 = single(filtfilt(sosF, gain, double(dataV6raw)'));
loadedV6chs = chnsNormal_V6(chp_V6);

%% plotting

figure('Name',figName, 'pos',[10 10 430 800]) %
subplot(2,1,2)

% define scaling (to micovolts)
scale2mV = 1/3276700*1E6;

plGap = 700*scale2mV; % this is the plotting gap between channel traces
skipN = 2; % samples to skip for more efficient plotting
plcolor = 'k';
lnwidth = 0.5;

xlimZoom = [428.1492 428.1492+0.2]; % 200ms zoom in
ylimZoom = [-118. 70]; % manual to have the same and fit both

ylimMult = [-2.48e+03 76.49]; % manual to have the same and fit both
xlimMult = [startRawT endRawT]; 

plot(timePlotRaw(1:skipN:end), scale2mV*dataV6(1:skipN:end,:)-((1:size(dataV6,2))-1)*plGap,'color', plcolor, 'linewidth', lnwidth)
hold on

box off
set(gca,'Color','w')
set(gca,'TickDir','out');
axis tight
xlabel('Time in session (s)')
set(gca,'ytick',[])
ax = gca;
ax.XColor = 'k';
ax.YColor = 'w';
axis off
ylim(ylimMult)

hold on
plot([xlimMult(end) xlimMult(end)], [0 200]+ylimMult(1), 'linewidth', 1, 'color', 'k')
text(xlimMult(end)+0.02,  100+ylimMult(1),'200 \muV','fontsize',5)

plot([xlimMult(end)-0.1 xlimMult(end)], [ylimMult(1) ylimMult(1)], 'linewidth', 1, 'color', 'k')
text(xlimMult(end)-0.06, ylimMult(1)-100,'0.1 s','fontsize',5)

% plot laser signal
h_spLED = subplot(2,1,1);
plot(TimeContLEDV(startLEDVI:endLEDVI), LEDV(startLEDVI:endLEDVI,1),'color',[0 1 0.8]*0.9, 'linewidth', 1.0)

axis off
axis tight
h_spLED.Position = [0.1300    0.465    0.7750    0.015];
set(gcf,'Color','w')
