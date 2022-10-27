
% Code to produce Supplementary Fig. 7
% Analysis of saccade amplitude and direcction during the optogenetic and visual detection taskS
% sessionIDs = [456:460]; % DC laser sessions with consistent stimulation

addpath sacc
pxPerDeg = 27.854; % pixels per degree of visual angle at the center


%% run saccade dection for all trials
% [sacc] = detectSaccEp(EYESeSc);
% save('SaccadeAnalysisData','sessionIDs', 'EYESeSc', 'PHOTeS', 'tr', 'sacc', 'pxPerDeg')

% alternatively, load data from file
load  Suppl_Fig7_data.mat


%% create labels for conditions and trial errors

trErrlabel{1} =  'Correct';
trErrlabel{2} =  'Maintained fix.';
trErrlabel{3} =  'False alarm';
trErrlabel{4} =  'Premature saccade';
trErrlabel{5} =  'No response to target';
trErrlabel{6} =  'Break fixation on target';
trErrlabel{7} =  'Break fixation';
trErrlabel{8} =  'Not facing screen';
trErrlabel{9} =  'No response to fix spot';

condlabel{1} =  'Catch';
condlabel{2} =  'High contr. visual';
condlabel{3} =  'High contr. visual + opto';
condlabel{4} =  'Low contr. visual';
condlabel{5} =  'Low contr. visual + opto';
condlabel{6} =  'Opto only';
condlabel{7} =  'Sham';

%% loop through every trial and try to extract saccade start and end points
preT = 0.6; % pre stimuglus onset analysis time
postT = 0.65+0.6; % post stimulus onset analysis time
timeXa = linspace(-preT,postT,size(EYESeSc,2)); % trial time
nTrls = numel(tr.trialNr);
thisRTsaccVal = nan(1,nTrls);
withinTrFAT = nan(1,nTrls);

for iT = 1:nTrls
    
    if mod(iT,100)==1
        disp(['Trial ' num2str(iT) '/' num2str(nTrls)])
    end
    
    thisTrError = tr.trialError(iT);
    thisTrCond = tr.condition(iT);
    
    if thisTrCond>1 && thisTrError == 1% correct go trials, non-catch
        fabydrift(iT) = false;
        withinTrRT = tr.reactionTime(iT)+tr.targetOn(iT)-tr.stimOnTphot(iT); % reaction time relative to 0 (stimOnT)
        isValidSacc = ~isnan(sacc.Offsets{iT}) & sacc.Offsets{iT}<=numel(timeXa) & ~isnan(sacc.Onsets{iT}) & sacc.Onsets{iT}<=numel(timeXa);
        sacOffT = timeXa(sacc.Offsets{iT}(isValidSacc));
        sacOnT = timeXa(sacc.Onsets{iT}(isValidSacc));
        
        if isempty(sacOffT) % no saccades detected, correct trial was caused by drift
            cobydrift(iT) = true;
            fasc.x(iT) = NaN;
            fasc.y(iT) = NaN;
            rtsc.x(iT) = NaN;
            rtsc.y(iT) = NaN;
        else
            cobydrift(iT) = false;
            
            % define saccade time as mean between on and offset
            saccMidT = mean([sacOffT; sacOnT]);
            
            % find the saccMidT closest to online RT (RT = time when fix window was left)
            [thisRTsaccVal(iT), thisRTsaccInd] = min(abs(saccMidT-withinTrRT));
            
            thisRTsaccOffT = sacOffT(thisRTsaccInd);
            thisRTsaccOnT = sacOnT(thisRTsaccInd);
            fromRTsaccOffI = find(timeXa >= thisRTsaccOffT,1,'first');
            fromRTsaccPreI = find(timeXa >= thisRTsaccOnT,1,'first');
            
            % take the median 25ms after the offset for x and y eye data
            rtsc.x(iT) = median(EYESeSc(iT,fromRTsaccOffI:fromRTsaccOffI+25,1));
            rtsc.y(iT) = median(EYESeSc(iT,fromRTsaccOffI:fromRTsaccOffI+25,2));
            rtsc.x_pre(iT) = median(EYESeSc(iT,fromRTsaccPreI-25:fromRTsaccPreI,1));
            rtsc.y_pre(iT) = median(EYESeSc(iT,fromRTsaccPreI-25:fromRTsaccPreI,2));
            fasc.x(iT) = NaN;
            fasc.y(iT) = NaN;
            fasc.x_pre(iT) = NaN;
            fasc.y_pre(iT) = NaN;
            
            % check if the detected saccade was away from the fix point
            saccFixDist(iT) = pdist([rtsc.x_pre(iT) rtsc.y_pre(iT); 0 0]) -  pdist([rtsc.x(iT) rtsc.y(iT); 0 0]);
            saccIsDirAway(iT) = saccFixDist(iT)<0;
            fabydrift(iT) = false;
            
            % if the closest saccade was more than 25ms away from the RT or if it was directed towards the fixation point, it was most likely caused by drift
            if thisRTsaccVal(iT)>0.025 || ~saccIsDirAway(iT)
                cobydrift(iT) = true;
                disp('Correct trial caused by drift.')
                rtsc.x(iT) = NaN;
                rtsc.y(iT) = NaN;
                fasc.x(iT) = NaN;
                fasc.y(iT) = NaN;
            end
        end
        
    elseif thisTrCond==1 && thisTrError ~= 1% incorrect catch trials = false alarms
        
        cobydrift(iT) = false;
        withinTrFAT(iT) = tr.falseAlarmTime(iT)+tr.acquireFixT(iT)-tr.stimOnTphot(iT); % reaction time relative to 0 (stimOnT)
        isValidSacc = ~isnan(sacc.Offsets{iT}) & sacc.Offsets{iT}<=numel(timeXa) & ~isnan(sacc.Onsets{iT}) & sacc.Onsets{iT}<=numel(timeXa);
        sacOffT = timeXa(sacc.Offsets{iT}(isValidSacc));
        sacOnT = timeXa(sacc.Onsets{iT}(isValidSacc));
        
        if isempty(sacOffT) % no saccades detected, false alam was caused by drift
            fabydrift(iT) = true;
            rtsc.x(iT) = NaN;
            rtsc.y(iT) = NaN;
            fasc.x(iT) = NaN;
            fasc.y(iT) = NaN;
        else
            
            fabydrift(iT) = false;
            
            % define saccade time as mean between on and offset
            saccMidT = mean([sacOffT; sacOnT]);
            
            % find the saccMidT closest to online RT (RT = time when fix window was left)
            [thisRTsaccVal(iT), thisRTsaccInd] = min(abs(saccMidT-withinTrFAT(iT)));
            
            thisRTsaccOffT = sacOffT(thisRTsaccInd);
            thisRTsaccOnT = sacOnT(thisRTsaccInd);
            fromRTsaccOffI = find(timeXa >= thisRTsaccOffT,1,'first');
            fromRTsaccPreI = find(timeXa >= thisRTsaccOnT,1,'first');
            
            % take the median 25ms after the offset for x and y eye data
            fasc.x(iT) = median(EYESeSc(iT,fromRTsaccOffI:fromRTsaccOffI+25,1));
            fasc.y(iT) = median(EYESeSc(iT,fromRTsaccOffI:fromRTsaccOffI+25,2));
            fasc.x_pre(iT) = median(EYESeSc(iT,fromRTsaccPreI-25:fromRTsaccPreI,1));
            fasc.y_pre(iT) = median(EYESeSc(iT,fromRTsaccPreI-25:fromRTsaccPreI,2));
            rtsc.x(iT) = NaN;
            rtsc.y(iT) = NaN;
            rtsc.x_pre(iT) = NaN;
            rtsc.y_pre(iT) = NaN;
            
            % check if the detected saccade was away from the fix point
            saccFixDist(iT) = pdist([fasc.x_pre(iT) fasc.y_pre(iT); 0 0]) -  pdist([fasc.x(iT) fasc.y(iT); 0 0]);
            saccIsDirAway(iT) = saccFixDist(iT)<0;
            cobydrift(iT) = false;
            
            % if the closest saccade was more than 25ms away from the RT or if it was directed towards the fixation point, it was most likely caused by drift
            if thisRTsaccVal(iT)>0.025 || ~saccIsDirAway(iT)
                fabydrift(iT) = true;
                disp('False alarm caused by drift.')
                rtsc.x(iT) = NaN;
                rtsc.y(iT) = NaN;
                fasc.x(iT) = NaN;
                fasc.y(iT) = NaN;
                
            end
        end
        
    else
        rtsc.x(iT) = NaN;
        rtsc.y(iT) = NaN;
        fasc.x(iT) = NaN;
        fasc.y(iT) = NaN;
        fabydrift(iT) = false;
        cobydrift(iT) = false;
    end
    
    
    [theta, rho] = cart2pol(fasc.x(iT),fasc.y(iT));
    
    if cobydrift(iT) || fabydrift(iT)
        disp('Drift or blink trial detected')
    end
    
    
    
end

nTrExluded = sum(fabydrift|cobydrift);
disp([num2str(nTrExluded) ' trials excluded because no valid saccade could be detected'])

%% convert x and y values in amplitude and direction (polar coordinates)
[rtsc.theta,rtsc.rho] = cart2pol(rtsc.x,rtsc.y);
[fasc.theta,fasc.rho] = cart2pol(fasc.x,fasc.y);
withinTrRT = tr.reactionTime+tr.targetOn-tr.stimOnTphot;

% calculate the target position
[target_theta, target_rho] = cart2pol(tr.circTargetPositionX(1)/pxPerDeg, tr.circTargetPositionY(1)/pxPerDeg);

%% Plot the results as categorical scatterplot
figure('Position', [200 200 550 1000])

% false alarm selection
fasel = tr.condition==1 & tr.trialError ~= 1 & ~fabydrift;
cdsel = tr.condition~=1 & tr.trialError == 1 & ~cobydrift;

vals_rho =    [rtsc.rho(tr.condition==2 & cdsel) rtsc.rho(tr.condition==3 & cdsel) rtsc.rho(tr.condition==4 & cdsel) rtsc.rho(tr.condition==5 & cdsel) rtsc.rho(tr.condition==6 & cdsel) rtsc.rho(tr.condition==7 & cdsel) fasc.rho(fasel)];
vals_theta =  [rtsc.theta(tr.condition==2 & cdsel) rtsc.theta(tr.condition==3 & cdsel) rtsc.theta(tr.condition==4 & cdsel) rtsc.theta(tr.condition==5 & cdsel) rtsc.theta(tr.condition==6 & cdsel) rtsc.theta(tr.condition==7 & cdsel) fasc.theta(fasel)];
groups = [2*ones(size(rtsc.rho(tr.condition==2 & cdsel))) 3*ones(size(rtsc.rho(tr.condition==3 & cdsel))) 4*ones(size(rtsc.rho(tr.condition==4 & cdsel))) 5*ones(size(rtsc.rho(tr.condition==5 & cdsel)))  6*ones(size(rtsc.rho(tr.condition==6 & cdsel))) 7*ones(size(rtsc.rho(tr.condition==7 & cdsel)))  1*ones(size(fasc.rho(fasel)))];

% get number of saccades per condition
nSaccCond(1) = sum(fasel); % number of false alarm saccades during catch trials
for iCond = 2:7
    nSaccCond(iCond) = sum(tr.condition==iCond & cdsel); % number of hit saccades during non-catch trials
end

% prepare condition labels with n=trials
condlabelWithN = condlabel(unique(groups));
for iCond = 1:7
    condlabelWithN{iCond} = [condlabelWithN{iCond} ' (' num2str(nSaccCond(iCond)) ')'];
end

subplot(2,1,1)
hold on
xt = [1 2 3.5 4.5 6 7 8 ]; % x centers for plotting
catScatterplot(vals_rho', groups', 'Labels', condlabelWithN,'MarkerSize', 7, 'BoxColor', [156, 95, 95]/255, 'plotOrder', [  2 3 4 5 6 1 7 ], 'xcenters', xt);
xtickangle(45)
ylabel('Saccade amplitude (°)')
xxx1 = xlim;
yyy = ylim;
ylim([0 yyy(2)])
yline(target_rho, 'm-') % location of target
% for i=RFecc; yline(i, 'k-'); end % shank 5 RF locations
set(gca,'TickDir','out');

subplot(2,1,2)
hold on

catScatterplot(rad2deg(vals_theta)', groups', 'Labels', condlabelWithN,'MarkerSize', 7,'BoxColor', [156, 95, 95]/255, 'plotOrder', [  2 3 4 5 6 1 7 ], 'xcenters', xt);
xtickangle(45)
ylabel('Saccade direction (°)')
ylim([-180 180])
xxx2 = xlim;
yline(rad2deg(target_theta), 'm-')% location of target
% for i=RFpol; yline(i, 'k-'); end % shank 5 RF locations
set(gca,'TickDir','out');

set(gcf, 'color','w')

% Plot RF location and size from shank 5 (example opto MUA shank with strong response)
RFecc =  [   12.2085   10.8024   10.9979];
RFeccMean = mean(RFecc);

RFpol = [ -70.4402  -67.5600  -67.5038];
RFpolMean = mean(RFpol);

RFsize = [ 7.4611    7.5748    7.3849];
RFsizeMean = mean(RFsize);

% calculate width of RF as polar angle span
RFpolSizeDeg= rad2deg(asin(RFsizeMean/(2*RFeccMean))*2);

subplot(2,1,1)
% yline(RFeccMean-0.5*RFsizeMean, 'k:')
% yline(RFeccMean+0.5*RFsizeMean, 'k:')
% yline(RFeccMean, 'k')
% for i=RFecc; yline(i); end
ylim([0 15.5])
plot([xxx1(2) xxx1(2)],[RFeccMean-0.5*RFsizeMean RFeccMean+0.5*RFsizeMean],'k', 'LineWidth', 1)


subplot(2,1,2)
% yline(RFpolMean-0.5*RFpolSizeDeg, 'k:')
% yline(RFpolMean+0.5*RFpolSizeDeg, 'k:')
% yline(RFpolMean, 'k')
% % for i=RFpol; yline(i); end
plot([xxx2(2) xxx2(2)],[RFpolMean-0.5*RFpolSizeDeg RFpolMean+0.5*RFpolSizeDeg],'k', 'LineWidth', 1)

% save plot
savePlotFlag = true;
if savePlotFlag
    set(gcf, 'Renderer', 'painters')
    figFileName = 'SaccadeAmplitudeAndAngle.pdf';
    print(gcf,figFileName,'-dpdf')
end

%% pairwise statistical test for saccade amplitudes
% Wilcoxon rank sum test is equivalent to the Mann-Whitney U test

% 'High contrast visual' vs. 'Low contrast visual'
[pp_amp(1), hh_amp(1), stats] = ranksum( rtsc.rho(tr.condition==2 & cdsel), rtsc.rho(tr.condition==4 & cdsel));

% 'High contrast' vs. 'High contrast + opto'
[pp_amp(2), hh_amp(2), stats] = ranksum( rtsc.rho(tr.condition==2 & cdsel), rtsc.rho(tr.condition==3 & cdsel));

% 'Low contrast' vs. 'Low contrast + opto'
[pp_amp(3), hh_amp(3), stats] = ranksum( rtsc.rho(tr.condition==4 & cdsel), rtsc.rho(tr.condition==5 & cdsel));

% 'Opto only' vs. 'Catch'
[pp_amp(4), hh_amp(4), stats] = ranksum( rtsc.rho(tr.condition==6 & cdsel), fasc.rho(fasel));

% 'Sham' vs. 'Catch'
[pp_amp(5), hh_amp(5), stats] = ranksum( rtsc.rho(tr.condition==7 & cdsel), fasc.rho(fasel));

% FDR correction
[h_adj_amp, crit_p_amp, adj_ci_cvrg_amp, adj_p_amp] = fdr_bh(pp_amp);



%% pairwise statistical test for saccade direction (angle)
% Watson's U2 statistic is similar to the Mann-Whitney U test but applied to cirular data

nPerms = 1000;
rng('default') % For reproducibility

% 'High contrast visual' vs. 'Low contrast visual'
[pp_dir(1), hh_dir(1), stats] = watsons_U2_perm_test( rtsc.theta(tr.condition==2 & cdsel)', rtsc.theta(tr.condition==4 & cdsel)', nPerms);

% 'High contrast' vs. 'High contrast + opto'
[pp_dir(2), hh_dir(2), stats] = watsons_U2_perm_test( rtsc.theta(tr.condition==2 & cdsel)', rtsc.theta(tr.condition==3 & cdsel)', nPerms);

% 'Low contrast' vs. 'Low contrast + opto'
[pp_dir(3), hh_dir(3), stats] = watsons_U2_perm_test( rtsc.theta(tr.condition==4 & cdsel)', rtsc.theta(tr.condition==5 & cdsel)', nPerms);

% 'Opto only' vs. 'Catch'
[pp_dir(4), hh_dir(4), stats] = watsons_U2_perm_test( rtsc.theta(tr.condition==6 & cdsel)', fasc.theta(fasel)', nPerms);

% 'Sham' vs. 'Catch'
[pp_dir(5), hh_dir(5), stats] = watsons_U2_perm_test( rtsc.theta(tr.condition==7 & cdsel)', fasc.theta(fasel)', nPerms);

% FDR correction
[h_adj_dir, crit_p_dir, adj_ci_cvrg_dir, adj_p_dir] = fdr_bh(pp_dir);


adj_p_dir