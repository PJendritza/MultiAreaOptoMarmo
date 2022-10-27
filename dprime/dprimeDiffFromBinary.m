function [ dpDiff, cDiff ] = dprimeDiffFromBinary( isHit1 ,isHit2, isFalseAlarm1, isFalseAlarm2 )
%DPRIMEFROMBINARY Summary of this function goes here
%   Detailed explanation goes here


hitRate1 = sum(isHit1)/size(isHit1,1);
falsAlarmtRate1 = sum(isFalseAlarm1)/size(isFalseAlarm1,1);
[dp1, c1] = getdprime(hitRate1,falsAlarmtRate1);

hitRate2 = sum(isHit2)/size(isHit2,1);
falsAlarmtRate2 = sum(isFalseAlarm2)/size(isFalseAlarm2,1);
[dp2, c2] = getdprime(hitRate2,falsAlarmtRate2);

dpDiff = dp1-dp2;
cDiff = c1-c2;

end