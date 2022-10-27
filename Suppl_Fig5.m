% Code to produce Supplementary Fig. 5
% Opto SUA population analysis and depth plot

addpath data

% load data
load Suppl_Fig5_data.mat

% define statistical thresholds
rsTalpha  = 0.05; % alpha threshold for ranksum test
stdThresh = 3.00; % sdt threshold over baseline

%% extract trial average firing rate for each unit
rate_cont = NaN(1,numel(cluSel));
rate_blue = NaN(1,numel(cluSel));
rate_cont_bl = NaN(1,numel(cluSel));
rate_blue_bl = NaN(1,numel(cluSel));

% epoch the SUA data
epocType = 'laserOn';
clu = epochSpikeTimes_clu(clu, epocType);

for iClu = 1:nSelClu
    
    thisCluID = selCluIDs(iClu);  % current cluster ID
    thisCluInd = selCluInds(iClu);  % current cluster index in data, not cluster ID!
    
    tr = clu.trCell{thisCluInd};
    trSel   = ones(size(tr.trialError));
    
    % check which animal
    thisMonkeyLetter = tr.monkeyLetter{1};
    monkeyLetter_iClu{iClu} = tr.monkeyLetter{1};
    
    % find opto vs control trials dependent on which monkey
    if strcmp(thisMonkeyLetter, 'A')
        trSel_blue = tr.thisLaserCond==7001 & tr.thisLaserSham==7010 & trSel; % 505nm opto condition
        trSel_yell = tr.thisLaserCond==7002 & tr.thisLaserSham==7010 & trSel; % 594nm control condition
        trSel_noLa = tr.thisLaserCond==7000 & tr.thisLaserSham==7010 & trSel; % no laser control condition
        
        trSel_opto = trSel_blue; % opto condition: 505nm
        trSel_cont = trSel_yell | trSel_noLa;
        
        nTr_opto_Monkey_A = sum(trSel_opto); % number of trials
        nTr_cont_Monkey_A = sum(trSel_cont); % number of trials
        
        % define opto stim and baseline times for statistical analysis
        optoStimFrom =  0.0;
        optoStimTo   =  0.25;
        optoStimDur = optoStimTo-optoStimFrom;
        
        optoBlFrom   = -0.25;
        optoBlTo     =  0;
        optoBlDur = optoBlTo-optoBlFrom;
        
    elseif strcmp(thisMonkeyLetter, 'D')
        trSel_blue = tr.thisLaserCond==7001 & tr.thisLaserSham==7010 & trSel; % 505nm opto condition
        % trSel_yell = tr.thisLaserCond==7002 & tr.thisLaserSham==7010 & trSel; % 594nm control condition was not used in this session for monkey D
        trSel_noLa = tr.thisLaserCond==7000 & tr.thisLaserSham==7010 & trSel; % no laser control condition
        
        trSel_opto = trSel_blue;
        trSel_cont = trSel_noLa;
        
        nTr_opto_Monkey_D = sum(trSel_opto); % number of trials
        nTr_cont_Monkey_D = sum(trSel_cont); % number of trials
        
        % define opto stim and baseline times for statistical analysis
        optoStimFrom =  0.0;
        optoStimTo   =  0.20;
        optoStimDur = optoStimTo-optoStimFrom;
        
        optoBlFrom   = -0.20;
        optoBlTo     =  0;
        optoBlDur = optoBlTo-optoBlFrom;
        
    else
        error('Could not identify animal!')
    end
    
    % control trials
    spkCell = clu.spikeTimesEp(trSel_cont,thisCluInd);
    
    % get number of spikes after stim
    spkMat = cell2mat(spkCell);
    nTr_cont = size(spkCell,1);
    rate_cont(thisCluInd) = sum(spkMat>optoStimFrom & spkMat<optoStimTo)/nTr_cont/optoStimDur; % stimulation
    rate_cont_bl(thisCluInd) = sum(spkMat>optoBlFrom & spkMat<optoBlTo)/nTr_cont/optoStimDur; % baseline
    
    % get number of spikes for each trial for stats
    for itr = 1:nTr_cont
        tmp =  spkCell{itr};
        rate_cont_tr{iClu}(itr) = sum(tmp>optoStimFrom & tmp<optoStimDur)/optoStimDur;
        rate_cont_bl_tr{iClu}(itr) = sum(tmp>optoBlFrom & tmp<optoBlTo)/optoStimDur;
    end
    
    % blue laser trials
    spkCell = clu.spikeTimesEp(trSel_blue,thisCluInd);
    
    % get number of spikes after opto stim
    spkMat = cell2mat(spkCell);
    nTr_blue = size(spkCell,1);
    rate_blue(thisCluInd) = sum(spkMat>optoStimFrom & spkMat<optoStimDur)/nTr_blue/optoStimDur;
    rate_blue_bl(thisCluInd) = sum(spkMat>optoBlFrom & spkMat<optoBlTo)/nTr_blue/optoStimDur; % baseline
    
    % get number of spikes for each trial for stats
    for itr = 1:nTr_blue
        tmp =  spkCell{itr};
        rate_blue_tr{iClu}(itr) = sum(tmp>optoStimFrom & tmp<optoStimDur)/optoStimDur;
        rate_blue_bl_tr{iClu}(itr) = sum(tmp>optoBlFrom & tmp<optoBlTo)/optoStimDur;
    end
    
    % get sdt increase over baseline
    stdIncr_blue(thisCluInd) = mean(rate_blue_tr{iClu})/std(rate_blue_bl_tr{iClu});
    
    % calculate rate increase (ratio) between opto and control
    rate_incr_blueVsCont{iClu} = mean(rate_blue_tr{iClu})/mean(rate_cont_tr{iClu});
    
    % ranksum test per unit - opto vs. control
    [p_blueVsCont(iClu),h_blueVsCont(iClu),stats] = ranksum(rate_blue_tr{iClu},rate_cont_tr{iClu},'alpha', rsTalpha);
    
    % identify significantly modulated channels
    isModulated_blueVsCont(iClu) = h_blueVsCont(iClu)==1 && stdIncr_blue(thisCluInd)>=stdThresh;
    
end

rate_bluestimbl = [rate_blue; rate_blue_bl];

% display number of trials per condition
disp(['Monkey A - opto trials:    n = ' num2str(nTr_opto_Monkey_A)])
disp(['Monkey A - control trials: n = ' num2str(nTr_cont_Monkey_A)])

disp(['Monkey D - opto trials:    n = ' num2str(nTr_opto_Monkey_D)])
disp(['Monkey D - control trials: n = ' num2str(nTr_cont_Monkey_D)])

%% diagonal plot - Opto (505nm) vs. Control (594nm & no laster)
figure('Position',[ 500  500  780  500])
subplot(1,2,1)
cla

mrkAlpha = 0.7;
mrkSize = 40;

for iClu = 1:nSelClu
    thisCluID = selCluIDs(iClu);  % current cluster ID!
    thisCluInd = selCluInds(iClu);  % current cluster index in data, not cluster ID!
    
    tr = clu.trCell{thisCluInd};
    thisMonkeyLetter = monkeyLetter_iClu{iClu};
    
    axmin = 0.1;
    axmax = max(rate_bluestimbl(:))+5;
    
    hold on
    if iClu ==1
        plot([axmin axmax], [axmin axmax], 'color', [0 0 0 0.7], 'linewidth', 0.5)
    end
    
    if strcmp(thisMonkeyLetter, 'A')
        mrkType = '^';
    else
        mrkType = 'o';
    end
    
    if isModulated_blueVsCont(iClu)
        mrkCol = [0 0.9 0.72];
        mrkEdgCol = 'k';
    else
        mrkCol = 'w';
        mrkEdgCol = 'k';
    end
    
    % scatter plot opto vs. control
    scatter(rate_cont(thisCluInd), rate_blue(thisCluInd), mrkSize, 'MarkerEdgeAlpha', mrkAlpha, 'MarkerFaceAlpha', mrkAlpha, 'MarkerFaceColor', mrkCol,'MarkerEdgeColor', mrkEdgCol, 'Marker', mrkType)
    
    axis([axmin axmax axmin axmax])
    
    ylabel('Opto SUA (spks/s)')
    xlabel('Control SUA (spks/s)')
    
    set(gca,'TickDir','out');
    set(gcf,'Color','w')
    drawnow
    
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    axis square
end

n_blueMod = sum(isModulated_blueVsCont);
pct_blueMod = n_blueMod*100/nSelClu;
disp(['Significantly modulated units : ' num2str(pct_blueMod) ' %'])

%% find electrode depth of single units
for iClu = 1:nSelClu
    
    thisCluID = selCluIDs(iClu);  % current cluster ID!
    thisCluInd = selCluInds(iClu);  % current cluster index in data, not cluster ID!
    
    % median subtract waveform
    medWf  = median(wfm.waveFormsTrimMean(thisCluInd,:,:),3);
    wfTmp = squeeze(wfm.waveFormsTrimMean(thisCluInd,:,:))-medWf';
    
    % find index of largest amplitude channel
    [maxV, maxInd] = max(abs(max(wfTmp,[],2) - min(wfTmp,[],2)));
    
    % y-coordinate = depth of largest amplitude channel
    depthClu(iClu) = clu.cfg.ycoords(maxInd);
    
    %     % plot
    %     imagesc(wfTmp)
    %     hold on
    %     plot(26, maxInd, 'xr')
    %     hold off
    %     drawnow
    %     pause
end

%% scatter plot of SUA accros depth of electrodes

% get the coordinates of the electrodes
xVec = clu.cfg.xcoords(1:32);
yVec = clu.cfg.ycoords(1:32);
yChVec = 1:32;

subplot(1,2,2)
hold on

for iClu = 1:nSelClu
    thisCluID = selCluIDs(iClu);  % current cluster ID!
    thisCluInd = selCluInds(iClu);  % current cluster index in data, not cluster ID!
    
    thisMonkeyLetter = monkeyLetter_iClu{iClu};
    
    if strcmp(thisMonkeyLetter, 'A')
        mrkType = '^';
    else
        mrkType = 'o';
    end
    
    if isModulated_blueVsCont(iClu)
        mrkCol = [0 0.9 0.72];
        mrkEdgCol = 'k';
    else
        mrkCol = 'w';
        mrkEdgCol = 'k';
    end
    
    scatter(rate_incr_blueVsCont{iClu}, depthClu(iClu), mrkSize, 'MarkerEdgeAlpha', mrkAlpha, 'MarkerFaceAlpha', mrkAlpha, 'MarkerFaceColor', mrkCol,'MarkerEdgeColor', mrkEdgCol, 'Marker', mrkType)
    
end

ylim([-850 50])
ylabel('Electrode depth (\mum)')
xlabel('^{Opto}/_{Control} SUA', 'Interpreter','tex')

set(gca,'TickDir','out');
set(gcf,'Color','w')
box off

set(gca, 'XScale', 'log')
xlim([0.6 20])
axis square
