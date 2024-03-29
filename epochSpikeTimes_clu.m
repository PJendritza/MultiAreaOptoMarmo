function [clu] = epochSpikeTimes_v2(clu, epocType)
% EPOCHSPIKETIMES Epoching of spike times into trials (=epochs) based on events from tr struct  
% Change to original version: take tr from clu.trCell per unit
%   Inputs:
%       clu.spikeTimes: {1×n cell}        % spike times in the session for each cluster 
%       tr.'field':     [1×trials double] % tr struct from loadSessions.m 
%       epocType:       % string for field name in tr e.g. 'stimOnTphot';    
%
%   Output:
%      clu.spikeTimesEp: {trials×n cell} % epoched spike times for each cluster and trial


nClusters = numel(clu.spikeTimes);

for u = 1:nClusters
    tr = clu.trCell{u};
    nTrials = numel(tr.trialNr);
    for t = 1:nTrials
        sel = clu.spikeTimes{1,u}>=tr.TrialStartT(t) & clu.spikeTimes{1,u} < tr.TrialEndT(t);
        clu.spikeTimesEp{t,u} = clu.spikeTimes{1,u}(sel)-tr.(epocType)(t); % shifted to event onset
%       clu.spikeTimesEp{t,u}(isnan(clu.spikeTimesEp{t,u}));
        clu.spikeTimesRawEp{t,u} = clu.spikeTimes{1,u}(sel); % not shifted
    end
end


end