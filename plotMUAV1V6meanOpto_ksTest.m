function [stats] = plotMUAV1V6meanOpto_ksTest( fullSessionName, tr, aMU1eS, aMU6eS, preT, postT, xlimPl, Fs, msSmooth)
% plots mean MUA psth across channels and makes ks test to find significantly modulated
% channels

sessionNr  = 475;
mapnames = {'V1_1A'; 'V1_1B'; 'V6_1A'; 'V6_1B'; 'V6_2A'; 'V6_2B'};

%% find sham trials
tr.isShamTr = tr.thisLaserSham ~= 7010;

%% identify responsive channels
timeXa = linspace(-preT,postT,size(aMU1eS,2)); % time for plotting
xSEM = 1; % factor to multiply the SEM for plotting

% define baseline and stimulus time windows
blstart   =  -0.25;
blend     =   0.00;
signstart =   0.00;
signsend  =   0.25;
stX = xlimPl(1)*1.05;  % x position for the plotting the asterisk

% trial average V1
MUAaTrAvrgV1 = squeeze(mean(aMU1eS));
% % baseline subtraction V1
MUAabaslineV1 = squeeze(mean(MUAaTrAvrgV1(find(timeXa>blstart,1,'first'):find(timeXa<blend,1,'last'),:)));
MUAaTrAvrgV1 = MUAaTrAvrgV1-MUAabaslineV1; % subtract mean
MUAabaslineV1std = squeeze(std(MUAaTrAvrgV1(find(timeXa>blstart,1,'first'):find(timeXa<blend,1,'last'),:)));
MUAaTrAvrgV1 = MUAaTrAvrgV1./MUAabaslineV1std; % devide by standare deviation

% % smooth  V1
MUAaTrAvrgV1 = gaussSmooth( MUAaTrAvrgV1', Fs, msSmooth);

% trial average V6
MUAaTrAvrgV6 = squeeze(mean(aMU6eS));
% % baseline subtraction V6
MUAabaslineV6 = squeeze(mean(MUAaTrAvrgV6(find(timeXa>blstart,1,'first'):find(timeXa<blend,1,'last'),:)));
MUAaTrAvrgV6 = MUAaTrAvrgV6-MUAabaslineV6; % subtract mean
MUAabaslineV6std = squeeze(std(MUAaTrAvrgV6(find(timeXa>blstart,1,'first'):find(timeXa<blend,1,'last'),:)));
MUAaTrAvrgV6 = MUAaTrAvrgV6./MUAabaslineV6std; % devide by standare deviation

% % smooth V6
MUAaTrAvrgV6 = gaussSmooth( MUAaTrAvrgV6', Fs, msSmooth);


% trial average V6 (505nm only)
MUAaTrAvrgV6_505nm = squeeze(mean(aMU6eS(tr.thisLaserCond == 7001 & ~tr.isShamTr, :,:)));
% % baseline subtraction V6
MUAabaslineV6_505nm = squeeze(mean(MUAaTrAvrgV6_505nm(find(timeXa>blstart,1,'first'):find(timeXa<blend,1,'last'),:)));
MUAaTrAvrgV6_505nm = MUAaTrAvrgV6_505nm-MUAabaslineV6_505nm; % subtract mean
MUAabaslineV6std_505nm = squeeze(std(MUAaTrAvrgV6_505nm(find(timeXa>blstart,1,'first'):find(timeXa<blend,1,'last'),:)));
MUAaTrAvrgV6_505nm = MUAaTrAvrgV6_505nm./MUAabaslineV6std_505nm; % devide by standare deviation

% % smooth V6
MUAaTrAvrgV6_505nm = gaussSmooth( MUAaTrAvrgV6_505nm', Fs, msSmooth);



%%% test if signal is significantly larger from baseline
ksAlpha = 0.05; % significance level for null hyp.
stdThresh = 3; % min. signal threshold

for cc = 1:size(MUAaTrAvrgV1,2)
    
    % define baseline and signal
    ksDataBaseli = aMU1eS(:, find(timeXa>blstart,1,'first'):find(timeXa<blend,1,'last'),cc);
    ksDataSignal = aMU1eS(:, find(timeXa>signstart,1,'first'):find(timeXa<signsend,1,'last'),cc);
    
    % find the maximum signal (in sigmas) in the stimulus period
    maxSign_V1(cc) = max(MUAaTrAvrgV1(find(timeXa>signstart,1,'first'):find(timeXa<signsend,1,'last'),cc));
    
    % test with ks statistic
    [h_V1(cc),  p_V1(cc)] = kstest2(ksDataBaseli(:), ksDataSignal(:),'Alpha',ksAlpha);
    
end

mod_V1 = h_V1 & abs(maxSign_V1)>stdThresh; % definition of "modulated channel"
disp([num2str(sum(mod_V1)) ' modulated V1 channels'])

for cc = 1:size(MUAaTrAvrgV6,2)
    
    % define baseline and signal
    ksDataBaseli = aMU6eS(:, find(timeXa>blstart,1,'first'):find(timeXa<blend,1,'last'),cc);
    ksDataSignal = aMU6eS(:, find(timeXa>signstart,1,'first'):find(timeXa<signsend,1,'last'),cc);
    
    % find the maximum signal (in sigmas) in the stimulus period
    maxSign_V6(cc) = max(MUAaTrAvrgV6(find(timeXa>signstart,1,'first'):find(timeXa<signsend,1,'last'),cc));
    
    % test with ks statistic
    [h_V6(cc),  p_V6(cc)] = kstest2(ksDataBaseli(:), ksDataSignal(:),'Alpha',ksAlpha);
    
end

mod_V6 = h_V6 & abs(maxSign_V6)>stdThresh; % definition of "modulated channel"
disp([num2str(sum(mod_V6)) ' modulated V6 channels'])

%% create trial and channel average plot

% pool all trials from modulated channels
tmp =  bsxfun(@minus, aMU1eS, reshape(MUAabaslineV1,1,1,[])); % subtract subtract mean baseline
tmp =  bsxfun(@rdivide, tmp, reshape(MUAabaslineV1std,1,1,[])); % devide by standard deviation

nConds = 2; % blue and yellow laser
trCondSel(1,:) = ~tr.isShamTr & tr.thisLaserCond == 7001 & tr.thisTargetAlpha ~= 0; % blue
trCondSel(2,:) = ~tr.isShamTr & tr.thisLaserCond == 7002 & tr.thisTargetAlpha ~= 0; % yellow

condColor(1,:) = [0.0 0.62 0.72]; % blue
condColor(2,:) = [230 148 0]/255; % yellow


for ishank = 1:2 % loop through the V1 shanks
    shankChns = zeros(1,64);
    shankChns(getMappingV1V6(mapnames{ishank},'normal', sessionNr)) =  1;
    for icond = 1:nConds
        % calculate mean across responsive channels
        MUAaTrPoolMeanV1{ishank,icond} = mean(reshape(permute(tmp(trCondSel(icond,:),:,mod_V1 & shankChns),[2,1,3]),size(aMU1eS,2),[],1),2);
        % calculate SEM across all pooled responsive channels and trials
        MUAaTrPoolStdV1{ishank,icond} = sem(reshape(permute(tmp(trCondSel(icond,:),:,mod_V1 & shankChns),[2,1,3]),size(aMU1eS,2),[],1)')';
        %calculate SEM for each responsive channel independenty
        for ich=find(mod_V1 & shankChns)
            ichBin = false(size(mod_V1)); ichBin(ich) = true;
            MUAaTrChSEMV1{ishank,icond,ich} = sem(reshape(permute(tmp(trCondSel(icond,:),:,ichBin & shankChns),[2,1,3]),size(aMU1eS,2),[],1)')';
        end
        % calculate mean SED across responsive channels
        try
            MUAaTrMeanSEMV1{ishank,icond} = mean(squeeze(cell2mat(MUAaTrChSEMV1(ishank,icond,mod_V1 & shankChns))),2);
        catch
            MUAaTrMeanSEMV1{ishank,icond} = [];
        end
    end
end


% pool all trials from modulated channels
tmp =  bsxfun(@minus, aMU6eS, reshape(MUAabaslineV6,1,1,[])); % subtract subtract mean baseline
tmp =  bsxfun(@rdivide, tmp, reshape(MUAabaslineV6std,1,1,[])); % devide by standard deviation

for ishank = 1:4  % loop through the V6 shanks
    shankChns = zeros(1,128);
    shankChns(getMappingV1V6(mapnames{ishank+2},'normal', sessionNr)) =  1;
    for icond = 1:nConds
        MUAaTrPoolMeanV6{ishank,icond} = mean(reshape(permute(tmp(trCondSel(icond,:),:,mod_V6 & shankChns),[2,1,3]),size(aMU6eS,2),[],1),2);
        MUAaTrPoolStdV6{ishank,icond} = sem(reshape(permute(tmp(trCondSel(icond,:),:,mod_V6 & shankChns),[2,1,3]),size(aMU6eS,2),[],1)')';
        %calculate SEM for each responsive channel independenty
        for ich=find(mod_V6 & shankChns)
            ichBin = false(size(mod_V6)); ichBin(ich) = true;
            MUAaTrChSEMV6{ishank,icond,ich} = sem(reshape(permute(tmp(trCondSel(icond,:),:,ichBin & shankChns),[2,1,3]),size(aMU6eS,2),[],1)')';
        end
        % calculate mean SED across responsive channels
        try
            MUAaTrMeanSEMV6{ishank,icond} = mean(squeeze(cell2mat(MUAaTrChSEMV6(ishank,icond,mod_V6 & shankChns))),2);
        catch
            MUAaTrMeanSEMV6{ishank,icond} = [];
        end
    end
end

%% plotting

figure('Name',fullSessionName, 'pos',[10 10 1600 800]) % V1 figure
timeXa = linspace(-preT,postT,size(MUAaTrAvrgV1,1));
set(gcf,'color','w');

subplot(4,4,9)
pmap = getMappingV1V6('V1_1A','normal', sessionNr);
imagesc(timeXa, 1:32, MUAaTrAvrgV1(:,pmap)');
colormap hot
hold on
plot([0 0], [1 32], '--w' )
xlabel('Time (s)')
% ylabel('Depth');
set(gca,'ytick',[])
title('V1 shank1A')
xlim([xlimPl])
set(gca,'TickDir','out');
% draw labels for significance testing
for i = find(mod_V1(pmap))
    text(stX, i,'*', 'color', [0 0 0],'FontSize',4,  'VerticalAlignment', 'middle')
end
box off
spPos = get(gca, 'Position');  % [left, bottom, width, height]
cb = colorbar;
cb.Position = [0.01+spPos(1)+spPos(3) spPos(2) 0.007 spPos(4)*0.95];


subplot(4,4,13)
timeXa = linspace(-preT,postT,size(MUAaTrAvrgV1,1));
% try plot(timeXa,gaussSmooth(MUAaTrAvrgV1(:,intersect(pmap,find(mod_V1)))',Fs, msSmooth), 'linewidth',1); catch; end
try
    hold on
    shadedErrorBar(timeXa,gaussSmooth( MUAaTrPoolMeanV1{1,2}', Fs, msSmooth), gaussSmooth( MUAaTrMeanSEMV1{1,2}'*xSEM, Fs, msSmooth),'lineprops',{'color', condColor(2,:),'linewidth',1},'transparent',1);
    shadedErrorBar(timeXa,gaussSmooth( MUAaTrPoolMeanV1{1,1}', Fs, msSmooth), gaussSmooth( MUAaTrMeanSEMV1{1,1}'*xSEM, Fs, msSmooth),'lineprops',{'color', condColor(1,:),'linewidth',1},'transparent',1);
catch
end
axis tight
xlim(xlimPl)
ylim(ylim+[-1 +1])

xlabel('Time (s)')
ylabel('MUA (z-score)')
set(gca,'TickDir','out');
box off


subplot(4,4,10)
pmap = getMappingV1V6('V1_1B','normal', sessionNr);
imagesc(timeXa, 1:32, MUAaTrAvrgV1(:,pmap)');
colormap hot
hold on
plot([0 0], [1 32], '--w' )
xlabel('Time (s)')
title('V1 shank1B')
xlim([xlimPl])
set(gca,'TickDir','out');
set(gca,'ytick',[])
box off
spPos = get(gca, 'Position');  % [left, bottom, width, height]
cb = colorbar;
cb.Position = [0.01+spPos(1)+spPos(3) spPos(2) 0.007 spPos(4)*0.95];

% draw labels for significance testing
for i = find(mod_V1(pmap))
    text(stX, i,'*', 'color', [0 0 0],'FontSize',4, 'VerticalAlignment', 'middle')
end

subplot(4,4,14)
timeXa = linspace(-preT,postT,size(MUAaTrAvrgV1,1));
% try plot(timeXa,gaussSmooth(MUAaTrAvrgV1(:,intersect(pmap,find(mod_V1)))',Fs, msSmooth), 'linewidth',1); catch; end
try
    hold on
    shadedErrorBar(timeXa,gaussSmooth( MUAaTrPoolMeanV1{2,2}', Fs, msSmooth), gaussSmooth( MUAaTrMeanSEMV1{2,2}'*xSEM, Fs, msSmooth),'lineprops',{'color', condColor(2,:),'linewidth',1},'transparent',1);
    shadedErrorBar(timeXa,gaussSmooth( MUAaTrPoolMeanV1{2,1}', Fs, msSmooth), gaussSmooth( MUAaTrMeanSEMV1{2,1}'*xSEM, Fs, msSmooth),'lineprops',{'color', condColor(1,:),'linewidth',1},'transparent',1);
catch
end
axis tight
xlim(xlimPl)
ylim(ylim+[-1 +1])
xlabel('Time (s)')
ylabel('MUA (z-score)')
set(gca,'TickDir','out');
box off


% figure('Name',fullSessionName, 'pos',[200 100 1500 400]) % V6 figure
% timeXa = linspace(-preT,postT,size(MUAaTrAvrgV6,1));
% set(gcf,'color','w');

subplot(4,4,1)
pmap = getMappingV1V6('V6_1A','normal', sessionNr);
imagesc(timeXa, 1:32, MUAaTrAvrgV6_505nm(:,pmap)');
colormap hot
hold on
plot([0 0], [1 32], '--w' )
xlabel('Time (s)')
% ylabel('Depth');
set(gca,'ytick',[])
title('V6 shank1A')
xlim([xlimPl])
set(gca,'TickDir','out');
% draw labels for significance testing
for i = find(mod_V6(pmap))
    text(stX, i,'*', 'color', [0 0 0],'FontSize',4, 'VerticalAlignment', 'middle')
end
box off
spPos = get(gca, 'Position');  % [left, bottom, width, height]
cb = colorbar;
cb.Position = [0.01+spPos(1)+spPos(3) spPos(2) 0.007 spPos(4)*0.95];

subplot(4,4,5)
timeXa = linspace(-preT,postT,size(MUAaTrAvrgV6_505nm,1));
% try plot(timeXa,gaussSmooth(MUAaTrAvrgV6(:,intersect(pmap,find(mod_V6)))',Fs, msSmooth), 'linewidth',1); catch; end
try
    hold on
    shadedErrorBar(timeXa,gaussSmooth( MUAaTrPoolMeanV6{1,2}', Fs, msSmooth), gaussSmooth( MUAaTrMeanSEMV6{1,2}'*xSEM, Fs, msSmooth),'lineprops',{'color', condColor(2,:),'linewidth',1},'transparent',1);
    shadedErrorBar(timeXa,gaussSmooth( MUAaTrPoolMeanV6{1,1}', Fs, msSmooth), gaussSmooth( MUAaTrMeanSEMV6{1,1}'*xSEM, Fs, msSmooth),'lineprops',{'color', condColor(1,:),'linewidth',1},'transparent',1);
catch
end
axis tight
xlim(xlimPl)
ylim(ylim+[-1 +1])
xlabel('Time (s)')
ylabel('MUA (z-score)')
set(gca,'TickDir','out');
box off

subplot(4,4,2)
pmap = getMappingV1V6('V6_1B','normal', sessionNr);
imagesc(timeXa, 1:32, MUAaTrAvrgV6_505nm(:,pmap)');
colormap hot
hold on
plot([0 0], [1 32], '--w' )
xlabel('Time (s)')
title('V6 shank1B')
xlim([xlimPl])
set(gca,'ytick',[])
set(gca,'TickDir','out');
% draw labels for significance testing
for i = find(mod_V6(pmap))
    text(stX, i,'*', 'color', [0 0 0],'FontSize',4, 'VerticalAlignment', 'middle')
end
box off
spPos = get(gca, 'Position');  % [left, bottom, width, height]
cb = colorbar;
cb.Position = [0.01+spPos(1)+spPos(3) spPos(2) 0.007 spPos(4)*0.95];

subplot(4,4,6)
timeXa = linspace(-preT,postT,size(MUAaTrAvrgV6_505nm,1));
% try plot(timeXa,gaussSmooth(MUAaTrAvrgV6(:,intersect(pmap,find(mod_V6)))',Fs, msSmooth), 'linewidth',1); catch; end
try
    hold on
    shadedErrorBar(timeXa,gaussSmooth( MUAaTrPoolMeanV6{2,2}', Fs, msSmooth), gaussSmooth( MUAaTrMeanSEMV6{2,2}'*xSEM, Fs, msSmooth),'lineprops',{'color', condColor(2,:),'linewidth',1},'transparent',1);
    shadedErrorBar(timeXa,gaussSmooth( MUAaTrPoolMeanV6{2,1}', Fs, msSmooth), gaussSmooth( MUAaTrMeanSEMV6{2,1}'*xSEM, Fs, msSmooth),'lineprops',{'color', condColor(1,:),'linewidth',1},'transparent',1);
catch
end
axis tight
xlim(xlimPl)
ylim(ylim+[-1 +1])
xlabel('Time (s)')
ylabel('MUA (z-score)')
set(gca,'TickDir','out');
box off

subplot(4,4,3)
pmap = getMappingV1V6('V6_2B','normal', sessionNr);
imagesc(timeXa, 1:32, MUAaTrAvrgV6_505nm(:,pmap)');
colormap hot
hold on
plot([0 0], [1 32], '--w' )
xlabel('Time (s)')
set(gca,'ytick',[])
title('V6 shank2B')
xlim([xlimPl])
set(gca,'TickDir','out');
% draw labels for significance testing
for i = find(mod_V6(pmap))
    text(stX, i,'*', 'color', [0 0 0],'FontSize',4, 'VerticalAlignment', 'middle')
end
box off
spPos = get(gca, 'Position');  % [left, bottom, width, height]
cb = colorbar;
cb.Position = [0.01+spPos(1)+spPos(3) spPos(2) 0.007 spPos(4)*0.95];

subplot(4,4,7)
timeXa = linspace(-preT,postT,size(MUAaTrAvrgV6_505nm,1));
% try plot(timeXa,gaussSmooth(MUAaTrAvrgV6(:,intersect(pmap,find(mod_V6)))',Fs, msSmooth), 'linewidth',1); catch; end
try
    hold on
    shadedErrorBar(timeXa,gaussSmooth( MUAaTrPoolMeanV6{4,2}', Fs, msSmooth), gaussSmooth( MUAaTrMeanSEMV6{4,2}'*xSEM, Fs, msSmooth),'lineprops',{'color', condColor(2,:),'linewidth',1},'transparent',1);
    shadedErrorBar(timeXa,gaussSmooth( MUAaTrPoolMeanV6{4,1}', Fs, msSmooth), gaussSmooth( MUAaTrMeanSEMV6{4,1}'*xSEM, Fs, msSmooth),'lineprops',{'color', condColor(1,:),'linewidth',1},'transparent',1);
catch
end
axis tight
xlim(xlimPl)
ylim(ylim+[-1 +1])
xlabel('Time (s)')
ylabel('MUA (z-score)')
set(gca,'TickDir','out');
box off

subplot(4,4,4)
pmap = getMappingV1V6('V6_2A','normal', sessionNr);
imagesc(timeXa, 1:32, MUAaTrAvrgV6_505nm(:,pmap)');
colormap hot
hold on
plot([0 0], [1 32], '--w' )
xlabel('Time (s)')
set(gca,'ytick',[])
title('V6 shank2A')
xlim([xlimPl])
set(gca,'TickDir','out');
% draw labels for significance testing
for i = find(mod_V6(pmap))
    text(stX, i,'*', 'color', [0 0 0],'FontSize',4, 'VerticalAlignment', 'middle')
end
box off
spPos = get(gca, 'Position');  % [left, bottom, width, height]
cb = colorbar;
cb.Position = [0.01+spPos(1)+spPos(3) spPos(2) 0.007 spPos(4)*0.95];

subplot(4,4,8)
timeXa = linspace(-preT,postT,size(MUAaTrAvrgV6_505nm,1));
% try plot(timeXa,gaussSmooth(MUAaTrAvrgV6(:,intersect(pmap,find(mod_V6)))',Fs, msSmooth), 'linewidth',1); catch; end
try
    hold on
    shadedErrorBar(timeXa,gaussSmooth( MUAaTrPoolMeanV6{3,2}', Fs, msSmooth), gaussSmooth( MUAaTrMeanSEMV6{3,2}'*xSEM, Fs, msSmooth),'lineprops',{'color', condColor(2,:),'linewidth',1},'transparent',1);
    shadedErrorBar(timeXa,gaussSmooth( MUAaTrPoolMeanV6{3,1}', Fs, msSmooth), gaussSmooth( MUAaTrMeanSEMV6{3,1}'*xSEM, Fs, msSmooth),'lineprops',{'color', condColor(1,:),'linewidth',1},'transparent',1);
catch
end
axis tight
xlim(xlimPl)
ylim(ylim+[-1 +1])
xlabel('Time (s)')
ylabel('MUA (z-score)')
set(gca,'TickDir','out');
box off



%% stats
% combine and save stats
stats.fullSessionName = fullSessionName;
stats.p_V1 = p_V1;
stats.p_V6 = p_V6;
stats.h_V1 = h_V1;
stats.h_V6 = h_V6;
stats.maxSign_V1 = maxSign_V1;
stats.maxSign_V6 = maxSign_V6;
stats.ksAlpha = ksAlpha;
stats.blstart   = blstart;
stats.blend     = blend;
stats.signstart = signstart;
stats.signsend  = signsend;
stats.msSmooth  = msSmooth;
stats.stdThresh = stdThresh;


end