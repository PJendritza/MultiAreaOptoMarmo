function [trlSTDSelection, trlChSTDfactor] = trlSTDSelec(aMU1eS, aMU6eS, cutoffSTD)
%trlSTDSelec
%   Code to find noisy trials in which the standard deviation of MUA across time was more than x-times larger than the median standard deviation across all epochs.

% Input:
% cutoffSTD     - desired threshold for standard deviation cutoff
% aMU1eS        - MUA in V1
% aMU6eS        - MUA in V6

stdTrials = squeeze(std(cat(3, aMU1eS, aMU6eS),1,2)); % get std of each trial to look for noisy trials

trlChSTDfactor = stdTrials./median(stdTrials); 
trlSTDSelection = logical(prod(trlChSTDfactor<cutoffSTD,2))'; % see if any of the channals crosses the threshold

disp([num2str(sum(~trlSTDSelection)) ' trials rejected'])
disp(['n = ' num2str(sum(trlSTDSelection)) ' trials'])


end

