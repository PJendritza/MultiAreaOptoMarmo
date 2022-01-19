% Fig7_fg
% Code to generate SUA examples to optogenetic stimulation and V6
% from Monkey D

% load data
addpath data
load('Fig7_fg_data','tr', 'cfg', 'clu','wfm', 'templates', 'winv', 'timeLEDe', 'ledSignal')

%% unwhiten all waveform templates and calculate amplitudes and depths
templatesUnW = zeros(size(templates));
for t = 1:size(templates,1)
    templatesUnW(t,:,:) = squeeze(templates(t,:,:))*winv;
end

% The amplitude on each channel is the positive peak minus the negative
clu.tempChanAmps = squeeze(max(templatesUnW,[],2))-squeeze(min(templatesUnW,[],2));

% The template amplitude is the amplitude of its largest channel (but see below for true tempAmps)
tempAmpsUnscaled = max(clu.tempChanAmps,[],2);

% need to zero-out the potentially-many low values on distant channels ...
threshVals = tempAmpsUnscaled*0.3; 
clu.tempChanAmps(bsxfun(@lt, clu.tempChanAmps, threshVals)) = 0;

% ... in order to compute the depth as a center of mass
cfg.templateDepths = sum(bsxfun(@times,clu.tempChanAmps,cfg.ycoords'),2)./sum(clu.tempChanAmps,2);


%% epoch cluster data
epocType = 'laserOn';
clu = epochSpikeTimes(clu,tr, epocType);

%% quality check for each clusters

% min number of spikes
clu.isLarge = clu.n_spikes>100; % find large clusters 
nClust = numel(clu.n_spikes);
nChs = size(cfg.chmapKS,1);

% percentage of spikes below the refractory period of 'cfg.msThres'
cfg.msThres = 1.5;
for i =1:nClust
clu.isiQuality(i) = sum(diff(clu.spikeTimes{i}*1000)<cfg.msThres)/numel(diff(clu.spikeTimes{i}*1000))*100; % in percent
end



%% Plotting part

% select the example single unit clusters from the list
selClu = cfg.cluIDList;

for iCluID =  selClu

figure('Name',cfg.fullSessionName, 'pos',[10 10 1200 800]) % new figure for each unit
set(gcf,'Color','w')

iCluInd = find(clu.ID==iCluID); % current cluster index in data, not cluster ID!

% get all spike times witing the trial (to exclude inter trial artifacts
tmp = clu.spikeTimesRawEp(:,iCluInd);
clu.spikeTimesRawEpV{1,iCluInd} = vertcat(tmp{:});

nChs = size(cfg.chmapKS,1);


%% autocorrelation

upFactor = 3; % upsample factor for autocorrelation

spkT_ms = round(clu.spikeTimes{iCluInd}*1000*upFactor); % in ms/x
tmpSignal = logical(zeros(1,max(spkT_ms)));
tmpSignal(spkT_ms) = true;

maxlags = 50*upFactor;
      [c,lags] = xcorr(tmpSignal,tmpSignal,maxlags);
      c(maxlags+1) = NaN; % remove trivial peak
      
      % normalize to max 0-1
       c = c/max(c);      

subplot(4,3,8); cla      
      plot(lags/upFactor,c, 'Color', [0 0 0], 'linewidth',0.5)
      box off

      xlabel('Lag (ms)')
      ylabel('Norm. autocorr.')
      set(gca,'TickDir','out');



      
%% waveform

% check which templates are used by the current cluster
thisTemplID = unique(clu.spikeTemplates(clu.spike_clusters == iCluID));
thisTemplInd = thisTemplID+1;

% check if the cluster only uses one template
if numel(thisTemplID)>1
    error('This cluster has more than one template.. not yet implemented')
end

% Find the amplitude and the largest channel for all templates
[tempAmpsUnscaled,  tempAmpMaxCh]= max(clu.tempChanAmps,[],2);

maxWfCh = tempAmpMaxCh(thisTemplInd); % this is the channel (KS channel mapping!) with the maximum amplitude in the template
maxWfCh_raw = cfg.chmapKS(maxWfCh)+1; % actual mapping from raw data, +1 for zero indexing of in ks

timeWf = ([1:wfm.nsamples])/25; % time in ms
timeWf = timeWf-mean(timeWf); % shift to center


%% all waveforms on relevant shank
h_spwf = subplot(3,3,[1,4,7]);
cla

allMeamWfs = squeeze(wfm.waveFormsTrimMean(iCluInd,:,:));

% subtract median to correct for offsets
allMeamWfs_n = allMeamWfs-median(allMeamWfs,2);
clu.allMeamWfs_n{iCluInd} = allMeamWfs_n; % copy into clu struct for saving

% find  amplitude
allMeamWfs_ampl = max(allMeamWfs_n,[],2)-min(allMeamWfs_n,[],2);

% find amplitude ranking
[~, allMeamWfs_amplRank] = sort(allMeamWfs_ampl); % max in ranking is max channel number

% find the highest amplitude channels and check how they spread in space
topChns = allMeamWfs_amplRank(end-1:end);
topChns_raw = cfg.chmapKS(topChns)+1;
% find their positions
xPosTop = cfg.xcoords(topChns);
yPosTop = cfg.ycoords(topChns);
% measure the dispersion
topSpread(iCluInd) = sqrt(std(xPosTop)^2 + std(yPosTop)^2);

if topSpread(iCluInd)<30
txtColor = 'k';   
else
txtColor = 'k';     
end

% scale to largest amplitude
allMeamWfs_nsc = allMeamWfs_n/max(allMeamWfs_ampl);

% find the shank of the unit

cla
hold on
for i = 1:size(allMeamWfs,1)
    chID = cfg.chmapKS(i)+1;
    isShankmember(i) = strcmp(cfg.shankMap.shankName{maxWfCh_raw} , cfg.shankMap.shankName{chID});
    if isShankmember(i)
        plot(110*timeWf+cfg.xcoords(i)*0, 50*(allMeamWfs_nsc(i,:))+cfg.ycoords(i),'color',[0 0 0], 'linewidth',0.5)
        % text(cfg.xcoords(i)+60, cfg.ycoords(i)-10, num2str(cfg.chmapKS(i)+1),'FontSize', 5, 'HorizontalAlignment', 'center', 'Color', txtColor)
    end
end

axis off
axis tight
% xlim([-550 550]-400)

h_spwf.Position = [0.32 0.1100 0.035 0.8150];

%% info text

subplot(4,3,3)
cla
hold off
axis off

% plot information about cluster as text in the same subplot
ofs = +10;
spc = -35;
text(0, ofs+10.5 , [cfg.base, ' - ' cfg.sessionType], 'FontSize', 10, 'Interpreter', 'none')
text(0, ofs+spc*1 ,['Cluster ID  =  ' num2str(iCluID)], 'FontSize', 10)
text(0, ofs+spc*2, ['n Spikes  =  ' num2str(clu.n_spikes(iCluInd))], 'FontSize', 10)
text(0, ofs+spc*3, ['Spikes<' num2str(cfg.msThres) 'ms ISI  =  ' num2str(clu.isiQuality(iCluInd)) '%'], 'FontSize', 10)
text(0, ofs+spc*4, ['Shank location  =  ' cfg.shankMap.shankName{maxWfCh_raw}], 'FontSize', 10, 'Interpreter', 'none')
text(0, ofs+spc*5, ['Depth  =  ' num2str(cfg.templateDepths(thisTemplInd)) '\mum'], 'FontSize', 10)

ylim([-200 10])

%% raster plot to stimulus
subplot(4,3,2) 
cla

if contains(epocType, 'laserOn')
    trSel_blue = tr.thisLaserCond==7001&tr.thisLaserSham==7010;
    trSel_yell = tr.thisLaserCond==7000&tr.thisLaserSham==7010;
    trSel_sham = tr.thisLaserCond==7001&tr.thisLaserSham==7011;
    trSel = trSel_blue;
    xlimPSTH = [-0.1 cfg.stimTime+0.1];
else
    trSel_blue = tr.thisLaserCond==7001;
    trSel_yell = ismember(tr.thisLaserCond, [7002 7000]);
    trSel_sham = tr.thisLaserCond==0;
    trSel = trSel_yell;
    xlimPSTH = [-0.2 1.0];
end

plotBlueLaserTrials = false; % manual flag to plot laser stim. trials

plotSpikeRaster( cellfun(@transpose,clu.spikeTimesEp(trSel,iCluInd),'UniformOutput',false),'PlotType','scatter');

set(gca, 'YDir','normal')
ylabel('Trial')
xlabel('Time (s)');
set(gcf,'Color','w')
set(gca,'TickDir','out');
h=get(gca);
h.XColor = 'none';
h.XAxis.Label.Color=[0 0 0];
h.XAxis.Label.Visible='on';
xlim(xlimPSTH)
box off


%%
if contains(epocType, 'laserOn')
    % draw laser
%     if ~exist('LEDVeS', 'var')
%         [~, ~, ~, ~, ~, ~, ~, ~, LEDVeS] =  loadSessions(cfg.sessionID, 'laserOn');
%     end
    preT = 0.6; % pre stimulus onset analysis time
    postT = 0.65+0.6; % post stimulus onset analysis time
    
%     timeLEDe = linspace(-preT, postT, size(LEDVeS,2));
%     ledSignal = mean(LEDVeS(trSel_blue,:,1))/max(mean(LEDVeS(trSel_blue,:,1)));
    hold on
    h_spLED = subplot(4,3,6);
    plot(timeLEDe, ledSignal, 'Color', [0 0.9 0.72])
    axis tight
    xlim(xlimPSTH)
    axis off
    h_spLED.Position = [0.4108    0.93    0.2134    0.015]; %h.Position; 
    
else
    % add stim indication
    yyy = ylim;
    hold on
    plot([0 cfg.stimTime], [yyy(2) yyy(2)]+5, 'k', 'LineWidth',1.5)
    ylim([yyy(1) yyy(2)+5])
end

%% PSTH plot to stimulus
subplot(4,3,5); hold on
cla

if contains(epocType, 'laserOn')
    msSmooth = 2; % smoothing window in ms
else
    msSmooth = 10; % smoothing window in ms
end

try
    % plot yellow laser, no sham
    spkCell = clu.spikeTimesEp(trSel_yell,iCluInd);
    % ph2 = plotSpikePSTH(spkCell,msSmooth, -0.5, 1.5);
    ph2 = SU_PSTH(spkCell,msSmooth, -0.5, 1.5);

end

if contains(epocType, 'laserOn')
    % plot blue laser, no sham
    spkCell = clu.spikeTimesEp(trSel_blue,iCluInd);
    % ph1 = plotSpikePSTH(spkCell,msSmooth, -0.5, 1.5);
    
    % set colors for laser psth
    ph2.mainLine.Color = [230 148 0]/255;
    ph2.edge(1).Color = ph2.mainLine.Color+(1-ph2.mainLine.Color)*0.55;
    ph2.edge(2).Color = ph2.mainLine.Color+(1-ph2.mainLine.Color)*0.55;
    ph2.patch.FaceColor = ph2.mainLine.Color;
    
    ph1 = SU_PSTH(spkCell,msSmooth, -0.5, 1.5);
    ph1.mainLine.Color = [0 0.9 0.72];
    ph1.edge(1).Color = ph1.mainLine.Color+(1-ph1.mainLine.Color)*0.55;
    ph1.edge(2).Color = ph1.mainLine.Color+(1-ph1.mainLine.Color)*0.55;
    ph1.patch.FaceColor = ph1.mainLine.Color;
    
    legend([ph1.mainLine; ph2.mainLine; ], [{'505 nm'}; {'495 nm'};],'location','northeast')
    legend boxoff
end

ylabel('Firing rate (spks/s)')
xlabel('Time (s)');
box off
xlim(xlimPSTH)
set(gca,'TickDir','out');


%% firing rate over session
subplot(4,3,11) 

spkCell = clu.spikeTimes(:,iCluInd); % {trials×1} cell array
spkCell{1} = spkCell{1}-cfg.startRawT; % normalize to start time
msSmooth = 2000; % smoothing window in ms

SU_PSTH(spkCell,msSmooth, 0, cfg.endRawT-cfg.startRawT);

box off
xlim([5 cfg.endRawT-cfg.startRawT-5])
xlabel('Recording time (s)')
ylabel('Firing rate (spks/s)')
% title('Session stability')
set(gca,'TickDir','out');

drawnow

    
end


