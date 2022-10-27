
function [CI_dp, CI_c] = getdprimeBootCI(nHits, nMiss, nFalseAlarms, nCorrectRejections)


nReps = 10000;
alpha = 0.05;        %alpha value

rng(1,'twister');

isHit        = [ones(1,nHits) zeros(1,nMiss)]';
isFalseAlarm = [ones(1,nFalseAlarms) zeros(1,nCorrectRejections)]';

n1 = numel(isHit);  %sample size 1
n2 = numel(isFalseAlarm); %sample size 2

    sampVectorX1 = isHit(ceil(rand(n1, nReps)*n1));
    sampVectorX2 = isFalseAlarm(ceil(rand(n2, nReps)*n2));
    [dpBoot, cBoot] = dprimeFromBinary(sampVectorX1,sampVectorX2);

CI_dp = prctile(dpBoot,[100*alpha/2,100*(1-alpha/2)]);
CI_c  = prctile(cBoot,[100*alpha/2,100*(1-alpha/2)]);



end


