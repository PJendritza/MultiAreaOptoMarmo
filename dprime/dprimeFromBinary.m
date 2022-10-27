function [ dp, c] = dprimeFromBinary( isHit , isFalseAlarm )
%DPRIMEFROMBINARY Calculate dprime and response bias from zeros and ones 
%   Detailed explanation goes here

hitRate = sum(isHit)/size(isHit,1);
falsAlarmtRate = sum(isFalseAlarm)/size(isFalseAlarm,1);

[dp, c] = getdprime(hitRate,falsAlarmtRate);


end