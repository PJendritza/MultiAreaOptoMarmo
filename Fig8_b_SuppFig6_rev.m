% Fig8_b
% Code to generate behavioral analysis for detection task with optogenetic stimulation

% load data
addpath data
load('Fig8_b_data','tr', 'sessionIDs')

%% recalculate the correct fals alarm rate
[ out, tr ] = recalcFalseAlarms( tr );

%% collect trial outcomes per condition
sessionSelection = sessionIDs;

conditions = [1:7];
clear hits misses falseAlarms breakFix invalidEyeResponse earlyEyeResponse initiatedTrials percentSacc nSacc pSacc binomialCI

for iCond = conditions
    selFreq = [ 0 ];
    
    hits{iCond}                 = sum(ismember(tr.sessionNr, sessionSelection) & ismember(tr.condition, iCond) & tr.flagStartT==1 & tr.trialError==1 & tr.gratcontrast~=0 & ismember(tr.thisLaserFreq, selFreq)  );
    misses{iCond}               = sum(ismember(tr.sessionNr, sessionSelection) & ismember(tr.condition, iCond) & tr.flagStartT==1 & tr.trialError==2 & tr.gratcontrast~=0 & ismember(tr.thisLaserFreq, selFreq) );
    falseAlarms{iCond}          = sum(ismember(tr.sessionNr, sessionSelection) & ismember(tr.condition, iCond) & tr.flagStartT==1 & tr.trialError==3 & tr.gratcontrast~=0 & ismember(tr.thisLaserFreq, selFreq) );
    breakFix{iCond}             = sum(ismember(tr.sessionNr, sessionSelection) & ismember(tr.condition, iCond) & tr.flagStartT==1 & tr.trialError==7  & ismember(tr.thisLaserFreq, selFreq) );
    invalidEyeResponse{iCond}	= sum(ismember(tr.sessionNr, sessionSelection) & ismember(tr.condition, iCond) & tr.flagStartT==1 & ismember(tr.trialError,[5 6])  & ismember(tr.thisLaserFreq, selFreq) );
    earlyEyeResponse{iCond}     = sum(ismember(tr.sessionNr, sessionSelection) & ismember(tr.condition, iCond) & tr.flagStartT==1 & tr.trialError==7  & ismember(tr.thisLaserFreq, selFreq) );
    initiatedTrials{iCond}      = sum(ismember(tr.sessionNr, sessionSelection) & ismember(tr.condition, iCond) & tr.flagStartT==1 & ismember(tr.trialError,[1:7])  & ismember(tr.thisLaserFreq, selFreq) );
    
    % stimulus onset time
    stOnT{iCond}                = tr.stimOnT(ismember(tr.sessionNr, sessionSelection) & ismember(tr.condition, iCond) & tr.flagStartT==1 & tr.trialError==3 & tr.gratcontrast~=0 & ismember(tr.thisLaserFreq, selFreq) );
    
    % reaction times
    rt{iCond}                   = tr.reactionTime(ismember(tr.sessionNr, sessionSelection) & ismember(tr.condition, iCond) & tr.flagStartT==1 & tr.trialError==1 & tr.gratcontrast~=0 & ismember(tr.thisLaserFreq, selFreq)  );
    rt_fa{iCond}                = tr.falseAlarmTime(ismember(tr.sessionNr, sessionSelection) & ismember(tr.condition, iCond) & tr.flagStartT==1 & tr.trialError==3 & tr.gratcontrast~=0 & ismember(tr.thisLaserFreq, selFreq) );
    
    
    
    if iCond==1
        percentSacc{1}  =  mean(out.falseAlarmRateRe) ;
        nSacc{1}  =  unique(out.nFalseAlarms +out.nCorrectReject) ;
        pSacc{1} = percentSacc{1}/100;
        binomialCI{1} = mean(out.binomialCI);
        
    else
        percentSacc{iCond}  = (hits{iCond}/(misses{iCond}+hits{iCond}))*100 ;
        nSacc{iCond}  =  hits{iCond} + misses{iCond} ;
        [pSacc{iCond}, binomialCI{iCond}] = binofit(hits{iCond},misses{iCond}+hits{iCond});
        
        
    end
    
    
    
end


% take the false alarm time minus stimulus onset time as reaction time for catch trials
rt{1} = rt_fa{1}; % CHECK THIS >>> !!!

barlabels={['Catch (' num2str(nSacc{1}) ')']; ['High constrast (' num2str(nSacc{2}) ')']; ['High constrast + opto (' num2str(nSacc{3}) ')']; ['Low constrast (' num2str(nSacc{4}) ')']; ['Low constrast + opto (' num2str(nSacc{5}) ')']; ['Opto only (' num2str(nSacc{6}) ')'] ; ['Sham (' num2str(nSacc{7}) ')']  };
strpX = [8 9];
flip6and7 = true;
% flip around condition 6 and 7 for a nicer plot
if flip6and7==true
    percentSacc = percentSacc([1 2 3 4 5 7 6]);
    binomialCI = binomialCI([1 2 3 4 5 7 6]);
    nSacc = nSacc([1 2 3 4 5 7 6]);
    barlabels = barlabels([1 2 3 4 5 7 6]);
    strpX = [8 9]-1;
end

figure('position',[500 50 300 500])

clf
xt = [1 2.5 3.5 5 6 7.5 8.5 ];
bh = bar(xt, cell2mat(percentSacc));
hold on


for pp = 1:100
    plot(strpX,[100 100-5]-pp*1.2,'w','linewidth',1)
end

for ic = 1:7
    % errorbar(xt(ic),pSacc{ic}*100, 100*(binomialCI{ic}(1)-pSacc{ic}), 100*(binomialCI{ic}(2)-pSacc{ic}),'k-')
    plot([xt(ic), xt(ic) ], [100*(binomialCI{ic}(1)), 100*(binomialCI{ic}(2)) ], 'k')
end

optoColor = [0 0.8 0.8]*0.9;
bh.FaceColor = 'flat';
bh.BarWidth = 0.9;
bh.CData(1,:) = [0.7 0.7 0.7];
bh.CData(2,:) = [0.4 0.4 0.4];
bh.CData(3,:) = optoColor;
bh.CData(4,:) = [0.4 0.4 0.4];
bh.CData(5,:) = optoColor;
bh.CData(6,:) = optoColor;
bh.CData(7,:) = optoColor;

bh.EdgeColor = 'none';


set(gca,'xticklabel',barlabels)
set(gca,'xTick',xt)
xtickangle(-60)
ylim([0 100])
ylabel('Saccade rate (%)')
% title('Behaviour - OptoGoNoGo task')
box off
set(gca,'TickDir','out');
set(gcf,'color','w');

ylim([23 100])


%% R code for pairwise comparisons of proportions with multiple comp. corr. (Benjamini Hochberg)
% see also Fig8_b_stats.R

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nhits  <- c(156.0520, 108.0000, 109.0000, 66.0000, 89.0000, 53.0000, 44.0000)
% ntrls <- c(467, 115, 123, 118, 109, 119, 116)
%
% pairwise.prop.test(nhits , ntrls, p.adjust.method = "BH")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run here: https://rdrr.io/snippets/

% Result:

% Pairwise comparisons using Pairwise comparison of proportions
%
% data:  nhits out of ntrls
%
%   1       2       3       4       5       6
% 2 < 2e-16 -       -       -       -       -
% 3 < 2e-16 0.250   -       -       -       -
% 4 1.9e-05 1.7e-10 5.7e-08 -       -       -
% 5 < 2e-16 0.013   0.221   9.5e-05 -       -
% 6 0.041   4.0e-15 2.6e-12 0.129   3.6e-08 -
% 7 0.419   < 2e-16 4.0e-15 0.013   1.7e-10 0.389
%
% P value adjustment method: BH


%% new reaction time plot

%% plot histogram of reaction times for each animal

figure('pos',[819   388   370   580])
clf
condlabels={['Catch (' num2str(nSacc{1}) ')']; ['High constrast (' num2str(nSacc{2}) ')']; ['High constrast + opto (' num2str(nSacc{3}) ')']; ['Low constrast (' num2str(nSacc{4}) ')']; ['Low constrast + opto (' num2str(nSacc{5}) ')']; ['Opto only (' num2str(nSacc{6}) ')'] ; ['Sham (' num2str(nSacc{7}) ')']  };

theColors(1,:) = [0.7 0.7 0.7];
theColors(2,:) = [0.4 0.4 0.4];
theColors(3,:) = [0 0.8 0.8]*0.9;
theColors(4,:) = [0.4 0.4 0.4];
theColors(5,:) = [0 0.8 0.8]*0.9;
theColors(6,:) = [0 0.8 0.8]*0.9;
theColors(7,:) = [0.4 0.4 0.4];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% high contrast
subplot(3,1,1)
for iCond = [2:3]
    
    
    hold on
    edges = [0:0.01:0.5];
    phh{iCond} = histogram(rt{iCond},edges,'EdgeColor', theColors(iCond,:), 'DisplayStyle','stairs', 'linewidth',0.5, 'linestyle','-');
    [f_rt{iCond}, rt_ksdens{iCond}] = ksdensity(rt{iCond});
    hold on
    phl{iCond} = plot(rt_ksdens{iCond}, f_rt{iCond}, 'Color', theColors(iCond,:),'linewidth',1.5, 'linestyle','-');
end
hold on
% histogram(RTsLowSF,80,'Normalization','probability','EdgeColor', 'none')
ylabel('Trials per bin')
xlabel('Reaction time (s)')
set(gcf,'color','white')
box off
xlim([0 0.5])
legend([phl{2}; phl{3}; ], [{'High contrast'}; {'High contrast + opto'};],'location','northeast')
legend boxoff
set(gca,'TickDir','out');

% Wilcoxon rank sum test is equivalent to the Mann-Whitney U test
[p_highContr, h_highContr, stats_highContr] = ranksum(rt{2},rt{3});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% low contrast
subplot(3,1,2)
for iCond = [4:5]
    hold on
    edges = [0:0.01:0.5];
    phh{iCond} = histogram(rt{iCond},edges,'EdgeColor', theColors(iCond,:), 'DisplayStyle','stairs', 'linewidth',0.5, 'linestyle','-');
    [f_rt{iCond}, rt_ksdens{iCond}] = ksdensity(rt{iCond});
    hold on
    phl{iCond} = plot(rt_ksdens{iCond}, f_rt{iCond}, 'Color', theColors(iCond,:),'linewidth',1.5, 'linestyle','-');
end
hold on
% histogram(RTsLowSF,80,'Normalization','probability','EdgeColor', 'none')
ylabel('Trials per bin')
xlabel('Reaction time (s)')
set(gcf,'color','white')
box off
xlim([0 0.5])
legend([phl{4}; phl{5}; ], [{'Low contrast'}; {'Low contrast + opto'};],'location','northeast')
legend boxoff
set(gca,'TickDir','out');

% Wilcoxon rank sum test is equivalent to the Mann-Whitney U test
[p_lowContr, h_lowContr, stats_lowContr] = ranksum(rt{4},rt{5});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% opto only vs sham
subplot(3,1,3)
for iCond = [6:7]
    hold on
    edges = [0:0.01:0.5];
    phh{iCond} = histogram(rt{iCond},edges,'EdgeColor', theColors(iCond,:), 'DisplayStyle','stairs', 'linewidth',0.5, 'linestyle','-');
    [f_rt{iCond}, rt_ksdens{iCond}] = ksdensity(rt{iCond});
    hold on
    phl{iCond} = plot(rt_ksdens{iCond}, f_rt{iCond}, 'Color', theColors(iCond,:),'linewidth',1.5, 'linestyle','-');
end
hold on
% histogram(RTsLowSF,80,'Normalization','probability','EdgeColor', 'none')
ylabel('Trials per bin')
xlabel('Reaction time (s)')
set(gcf,'color','white')
box off
xlim([0 0.5])

legend([phl{6}; phl{7}; ], [{'Opto only'}; {'Sham'};],'location','northeast')
legend boxoff
set(gca,'TickDir','out');

% Wilcoxon rank sum test is equivalent to the Mann-Whitney U test
[p_optoOnly, h_optoOnly, stats_optoOnly] = ranksum(rt{6},rt{7});


%% statistical results for reaction time analysis

ppVal = [p_highContr p_lowContr p_optoOnly]; % combine pvalues

% % FDR correction 
% [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(ppVal);

%% dprime analysis
addpath dprime


saccRate = cell2mat(percentSacc);
nHitsV = (saccRate/100).*cell2mat(nSacc);

% low contrast visual vs. catch
hitrate = saccRate(4)/100;
farate = saccRate(1)/100;
[dp_lowContVsCatch,c_lowContVsCatch] = getdprime(hitrate,farate);

% get bootstrap confidence intervals
nHits1 = hitrate*nSacc{4};
nMiss1 = round((1-hitrate)*nSacc{4});
nFalseAlarms1 = round(farate*nSacc{1}); % integer required for bootstrapping
nCorrectRejections1 = round((1-farate)*nSacc{1}); % integer required for bootstrapping
[CI_dp_lowContVsCatch, CI_c_lowContVsCatch] = getdprimeBootCI(nHits1, nMiss1, nFalseAlarms1, nCorrectRejections1);


% low contrast + opto vs. opto only
hitrate = saccRate(5)/100;
farate = saccRate(7)/100;
[dp_lowContOptoVsOptoOnly,c_lowContOptoVsOptoOnly] = getdprime(hitrate,farate);

% get bootstrap confidence intervals
nHits2 = hitrate*nSacc{5};
nMiss2 = round((1-hitrate)*nSacc{5});
nFalseAlarms2 = round(farate*nSacc{7}); % integer required for bootstrapping
nCorrectRejections2 = round((1-farate)*nSacc{7}); % integer required for bootstrapping
[CI_dp_lowContOptoVsOptoOnly, CI_c_lowContOptoVsOptoOnly] = getdprimeBootCI(nHits2, nMiss2, nFalseAlarms2, nCorrectRejections2);

figure
subplot(1,2,1); hold on
bar([1, 2],[dp_lowContVsCatch dp_lowContOptoVsOptoOnly])
plot([1, 1 ], [CI_dp_lowContVsCatch(1), CI_dp_lowContVsCatch(2)], 'k')
plot([2, 2 ], [CI_dp_lowContOptoVsOptoOnly(1), CI_dp_lowContOptoVsOptoOnly(2)], 'k')
set(gca, 'XTick', [1 2]);
set(gca,'xticklabel',{'Low contrast vs. Catch','Low contrast + opto vs. Opto only'})
ylabel('d''')
title({'Sensitivity';''})
set(gcf,'color','white')
set(gca,'TickDir','out');
box off
xtickangle(45)

subplot(1,2,2); hold on
bar([1, 2],[c_lowContVsCatch c_lowContOptoVsOptoOnly])
plot([1, 1 ], [CI_c_lowContVsCatch(1), CI_c_lowContVsCatch(2)], 'k')
plot([2, 2 ], [CI_c_lowContOptoVsOptoOnly(1), CI_c_lowContOptoVsOptoOnly(2)], 'k')
set(gca, 'XTick', [1 2]);
set(gca,'xticklabel',{'Low contrast vs. Catch','Low contrast + opto vs. Opto only'})
ylabel('c')
title({'Response bias';''})
set(gcf,'color','white')
set(gca,'TickDir','out');
box off
xtickangle(45)

% Compute bootstrap confidence intervals and p-value for the difference between the two d-prime (and response bias) values

[CI_dpDiff, p_dpDiff, CI_cDiff, p_cDiff] = getdprimeDiffBootCI(nHits2, nHits1, nMiss2, nMiss1, nFalseAlarms2, nFalseAlarms1, nCorrectRejections2 , nCorrectRejections1);

subplot(1,2,1)
yylim = ylim;
xxticks = xticks;
text(mean(xxticks), max(yylim)+0.04,['p = ' num2str(p_dpDiff,3)],'VerticalAlignment','bottom', 'HorizontalAlignment', 'center')
ylim([yylim(1) max(yylim)+0.1])

subplot(1,2,2)
yylim = ylim;
xxticks = xticks;
text(mean(xxticks), max(yylim)+0.04,['p <= ' num2str(p_cDiff,3)],'VerticalAlignment','bottom', 'HorizontalAlignment', 'center')
ylim([yylim(1) max(yylim)+0.1])