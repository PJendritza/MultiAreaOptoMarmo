
% Figure 5_a
% Code to generate spiking example traces from V1 and V6 data

%% laod raw data
addpath data
load('Fig5_a_data.mat')

% start end end time when data was recorded
sessionStartT = 298;
sessionEndT = 939;

%% define which channels to plot

% for V1 plot
chp_V1 = [14:25];   % plot these channels  from the shank
chnsNormal_V1 = [42,40,44,46,10,8,34,36,4,6,38,48,2,12,16,14,30,28,60,50,32,22,18,20,56,58,54,52,24,26,64,62]; % channel mapping of the shank

% for V6 plot
chp_V6 = [12:23]; % plot these channels  from the shank
chnsNormal_V6 = [61,59,29,23,63,53,49,51,17,31,19,21,55,57,25,27,15,1,13,11,41,39,7,5,35,37,3,9,33,43,47,45]; % channel mapping of the shank

%% define time vector
FsRaw = 25000; % raw sampling rate

% set start and end times manually
startRawT = 424;
endRawT = startRawT+10;

startRawI = (startRawT-sessionStartT)*FsRaw;
endRawI =   (endRawT-sessionStartT)*FsRaw-1;

timePlotRaw = linspace(startRawT,endRawT, endRawI-startRawI+1);

%% filter data
HighPass = 300;
LowPass = 6000;

[Z,P,K] = butter(4, [HighPass LowPass]/(FsRaw/2));
[sosF, gain] = zp2sos(Z,P,K);
disp('Filtering with 4th order Butterworth filter (300-6000Hz)')

dataV1 = single(filtfilt(sosF, gain, double(dataV1raw)'));
loadedV1chs = chnsNormal_V1(chp_V1);

dataV6 = single(filtfilt(sosF, gain, double(dataV6raw)'));
loadedV6chs = chnsNormal_V6(chp_V6);


%% prepare for plotting

% define scaling (to micovolts)
scale2mV = 1/3276700*1E6;

plGap = 700*scale2mV; % this is the plotting gap between channel traces

skipN = 2; % samples to skip for more efficient plotting

plcolor = 'k';
lnwidth = 0.5;

chZoomV1 = 28; % ch to zoom in as seen in plot
chZoomV6 = 11; % ch to zoom in as seen in plot

xlimZoom = [428.1492 428.1492+0.2]; % 200ms zoom in window
ylimZoom = [-118. 70]; % plotting ylimits

ylimMult = [-2.48e+03 76.49]; % plotting ylimits
xlimMult = [startRawT endRawT];


%% open figure and plot

figure('Name','Fig5a', 'pos',[10 10 800 500]) %
hold on

% V1 plot
p2 = subplot(3,2,3);
plot(timePlotRaw(1:skipN:end), scale2mV*dataV1(1:skipN:end,:)-((1:size(dataV1,2))-1)*plGap,'color', plcolor, 'linewidth', lnwidth)
hold on
yPos_V1 = -(find(loadedV1chs==chZoomV1)-1)*plGap;
plot(xlimZoom, [yPos_V1 yPos_V1], 'r.')

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

% V6 plot
p3 = subplot(3,2,4);
plot(timePlotRaw(1:skipN:end), scale2mV*dataV6(1:skipN:end,:)-((1:size(dataV6,2))-1)*plGap,'color', plcolor, 'linewidth', lnwidth)
hold on
yPos_V6 = -(find(loadedV6chs==chZoomV6)-1)*plGap;
plot(xlimZoom, [yPos_V6 yPos_V6], 'r.')

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
text(xlimMult(end)+0.1,  100+ylimMult(1),'200 \muV','fontsize',4)

plot([xlimMult(end)-0.5 xlimMult(end)], [ylimMult(1) ylimMult(1)], 'linewidth', 1, 'color', 'k')
text(xlimMult(end)-0.5, ylimMult(1)-100,'0.5 s','fontsize',4)

% V1 zoom in plot
p1 = subplot(3,2,1:2);
plot(timePlotRaw(:), scale2mV*dataV1(:,find(loadedV1chs==chZoomV1)),'color', plcolor, 'linewidth', lnwidth)
box off
set(gca,'Color','w')
set(gca,'TickDir','out');
axis tight
set(gca,'ytick',[])
ax = gca;
ax.XColor = 'k';
ax.YColor = 'w';
axis off
xlim(xlimZoom)
ylim(ylimZoom)

% V6 zoom in plot
p4 = subplot(3,2,5:6);
plot(timePlotRaw(:), scale2mV*dataV6(:,find(loadedV6chs==chZoomV6)),'color', plcolor, 'linewidth', lnwidth)
box off
set(gca,'Color','w')
set(gca,'TickDir','out');
axis tight
set(gca,'ytick',[])
ax = gca;
ax.XColor = 'k';
ax.YColor = 'w';
axis off
xlim(xlimZoom)
ylim(ylimZoom)

% add scale bars
hold on
plot([xlimZoom(end) xlimZoom(end)], [-50 0]-50, 'linewidth', 2, 'color', 'k')
text(xlimZoom(end)+0.002, -50-25,'50 \muV')

plot([xlimZoom(end)-0.01 xlimZoom(end)], [ylimZoom(1) ylimZoom(1)], 'linewidth', 2, 'color', 'k')
text(xlimZoom(end)-0.011, ylimZoom(1)-25,'10 ms')

% general plot attributes
set(gcf,'Color','w')

% set subplot positions
p1.Position = [0.1300    0.7093    0.7750    0.2157]; % [left bottom width height]
p2.Position = [0.1300    0.34    0.37    0.34];
p3.Position = [0.53    0.34    0.37    0.34];
p4.Position = [0.1300    0.1100    0.7750    0.2157];



