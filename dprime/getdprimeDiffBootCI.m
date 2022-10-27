function [CI_dpDiffBoot, p_dpDiffBoot, CI_cDiffBoot, p_cDiffBoot] = getdprimeDiffBootCI(nHits1,nHits2, nMiss1, nMiss2, nFalseAlarms1, nFalseAlarms2, nCorrectRejections1 , nCorrectRejections2)

nReps = 10000;
alpha = 0.05;        %alpha value

rng(1,'twister');

isHit1        = [ones(1,nHits1) zeros(1,nMiss1)]';
isFalseAlarm1 = [ones(1,nFalseAlarms1) zeros(1,nCorrectRejections1)]';
n1_1 = numel(isHit1);  %sample size 1
n2_1 = numel(isFalseAlarm1); %sample size 2

isHit2        = [ones(1,nHits2) zeros(1,nMiss2)]';
isFalseAlarm2 = [ones(1,nFalseAlarms2) zeros(1,nCorrectRejections2)]';
n1_2 = numel(isHit2);  %sample size 1
n2_2 = numel(isFalseAlarm2); %sample size 2

    sampVectorX1_1 = isHit1(ceil(rand(n1_1, nReps)*n1_1));
    sampVectorX2_1 = isFalseAlarm1(ceil(rand(n2_1, nReps)*n2_1));
    
    sampVectorX1_2 = isHit2(ceil(rand(n1_2, nReps)*n1_2));
    sampVectorX2_2 = isFalseAlarm2(ceil(rand(n2_2, nReps)*n2_2));
    
    % [ dpDiff ] = dprimeDiffFromBinary( isHit1 ,isHit2, isFalseAlarm1, isFalseAlarm2 )
    [dpDiffBoot  , cDiffBoot ] = dprimeDiffFromBinary(sampVectorX1_1, sampVectorX2_1, sampVectorX1_2, sampVectorX2_2);

CI_dpDiffBoot = prctile(dpDiffBoot,[100*alpha/2,100*(1-alpha/2)]);
CI_cDiffBoot = prctile(cDiffBoot,[100*alpha/2,100*(1-alpha/2)]);

% get p-values from bootstrap distributions
p_dpDiffBoot = sum((dpDiffBoot-mean(dpDiffBoot)) <= -abs(mean(dpDiffBoot)) | (dpDiffBoot-mean(dpDiffBoot)) >= abs(mean(dpDiffBoot))) / nReps;
if p_dpDiffBoot==0;  p_dpDiffBoot = 1/nReps; end % lower bound for p-value

p_cDiffBoot  = sum((cDiffBoot-mean(cDiffBoot)) <= -abs(mean(cDiffBoot)) | (cDiffBoot-mean(cDiffBoot)) >= abs(mean(cDiffBoot))) / nReps;
if p_cDiffBoot==0;  p_cDiffBoot = 1/nReps; end  % lower bound for p-value
end






