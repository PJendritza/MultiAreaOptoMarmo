function [y] = sem(x)
%SEM Calculation of the standard error of the mean
%   For vectors, Y = SEM(X) returns the standard error of the mean. For
%   N-D arrays, SEM operates along the first non-singleton dimension of X.

dim = find(size(x) ~= 1,1,'first'); % find first non-singleton dimension
y = std(x,0,dim)/sqrt(size(x,dim));

end

