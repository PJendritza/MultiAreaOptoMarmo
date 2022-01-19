function [ out, tr ] = recalcFalseAlarms( tr )

%% script to recalculate actual false alarm rate
% based on using the timing from all available Go trials and checking if
% for any catch trial the monkey would have done a false alarm or not

%% Calculate the time limits for every Go trial, defined after stimOn (t=0)

for iTr = 1:numel(tr.trialNr)
    
    if tr.condition(iTr)~=1
        tr.BRKstartT(iTr) = 0;
        tr.BRKendT(iTr)   = tr.stimTime(iTr);
        
        tr.PREstartT(iTr) = tr.stimTime(iTr);
        tr.PREendT(iTr)   = tr.stimTime(iTr)+50;
        
        tr.HITstartT(iTr) = tr.stimTime(iTr)+50;
        tr.HITendT(iTr)   = tr.stimTime(iTr)+500;
        
        tr.MISstartT(iTr) = tr.stimTime(iTr)+500;
        tr.MISendT(iTr)   = 1000*tr.TrialEndT(iTr)-1000*tr.TrialStartT(iTr)-tr.holdFix(iTr);
    else
        tr.BRKstartT(iTr) = NaN;
        tr.BRKendT(iTr)   = NaN;
        
        tr.PREstartT(iTr) = NaN;
        tr.PREendT(iTr)   = NaN;
        
        tr.HITstartT(iTr) = NaN;
        tr.HITendT(iTr)   = NaN;
        
        tr.MISstartT(iTr) = NaN;
        tr.MISendT(iTr)   = NaN;
    end
    
end


%% pick a random catch trial and random go trial
% Look at every catch (noGo) trial and see if it would have been a Hit or Miss

goTrNums = find(tr.condition~=1 & tr.flagStartT==1); % all go trials, which ones to take?
catchTrNums = find(tr.condition==1 & tr.flagStartT==1); % all go trials, which ones to take?

nHits = sum(tr.trialError==1 & tr.condition~=1 & tr.flagStartT==1);% total number of hits in dataset
nMiss = sum(tr.trialError==2 & tr.condition~=1 & tr.flagStartT==1);% total number of missed in dataset
nGoTr = nHits+nMiss; % total number of go trials in dataset

nCatchTrlsShuff = round(nGoTr/0.6*0.4); % total number of catch trials for shuffling: 40% of all trials should be catch trials (in task code)

out.nShuff = 1000; % number of shuffling iterations

rng(1); % set random number generator for reproducibility

for iShuff = 1:out.nShuff
    
    
    faRe.isBRK{iShuff}  = false(nCatchTrlsShuff,1);
    faRe.isPRE{iShuff}  = false(nCatchTrlsShuff,1);
    faRe.isHIT{iShuff}  = false(nCatchTrlsShuff,1);
    faRe.isMIS{iShuff}  = false(nCatchTrlsShuff,1);
    
    iTr = 1;
    totalNum = 0;
    while totalNum < nCatchTrlsShuff
        
        randInGo = randperm(length(goTrNums),1);   % generate random index for Go trial
        whichGo{iShuff} =  goTrNums(randInGo); % take a ramdom trial from all go trials
        
        randInCatch = randperm(length(catchTrNums),1);   % generate random index for Go trial
        whichCatch{iShuff} =  catchTrNums(randInCatch); % take a ramdom trial from all go trials
        
        
        if tr.falseAlarmTime(whichCatch{iShuff})*1000 < tr.HITstartT(whichGo{iShuff})
            faRe.isBRK{iShuff}(iTr) = true; % break fix
        end
        if tr.falseAlarmTime(whichCatch{iShuff})*1000 >= tr.PREstartT(whichGo{iShuff}) && tr.falseAlarmTime(whichCatch{iShuff})*1000 < tr.PREendT(whichGo{iShuff})
            faRe.isPRE{iShuff}(iTr) = true; % premature resp
        end
        if tr.falseAlarmTime(whichCatch{iShuff})*1000 >= tr.HITstartT(whichGo{iShuff}) && tr.falseAlarmTime(whichCatch{iShuff})*1000 < tr.HITendT(whichGo{iShuff})
            faRe.isHIT{iShuff}(iTr) = true; % hit
        end
        
        if tr.falseAlarmTime(whichCatch{iShuff})*1000 >= tr.MISstartT(whichGo{iShuff})
            faRe.isMIS{iShuff}(iTr) = true; % miss
        end
        
        
        totalNum = sum(faRe.isHIT{iShuff})+sum(faRe.isMIS{iShuff});
        iTr = iTr+1;
    end
    
    out.nFalseAlarms(iShuff) = sum(faRe.isHIT{iShuff});
    out.nCorrectReject(iShuff) = sum(faRe.isMIS{iShuff});
    
    out.falseAlarmRateRe(iShuff) = out.nFalseAlarms(iShuff)*100/(out.nFalseAlarms(iShuff)+out.nCorrectReject(iShuff));
    out.whichCatch = whichCatch;
    out.whichGo = whichGo;
    
    [out.pSacc(iShuff), out.binomialCI(iShuff,:)] = binofit(out.nFalseAlarms(iShuff), out.nFalseAlarms(iShuff) + out.nCorrectReject(iShuff));
    
    if mod(iShuff,100)==1
    disp(['Shuffling iteration: ' num2str(iShuff) '/' num2str(out.nShuff)])
    end
    
end


end