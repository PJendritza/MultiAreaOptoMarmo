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

barlabels={['Catch (' num2str(nSacc{1}) ')']; ['High constrast (' num2str(nSacc{2}) ')']; ['High constrast + opto (' num2str(nSacc{3}) ')']; ['Low constrast (' num2str(nSacc{4}) ')']; ['Low constrast + opto (' num2str(nSacc{5}) ')']; ['Opto only (' num2str(nSacc{6}) ')'] ; ['Sham (' num2str(nSacc{7}) ')']  };
strpX = [8 9];
flip6and7 = true;
% flip around condition 6 and 7 for a nicer plot
if flip6and7==true
    percentSacc = percentSacc([1 2 3 4 5 7 6]);
    binomialCI = binomialCI([1 2 3 4 5 7 6]);
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


