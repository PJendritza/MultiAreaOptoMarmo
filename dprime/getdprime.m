function [dp,c] = getdprime(h,fA)
% dprime d' calculation from hit and false alarm rate
%   outputs: 
%   dp = d'
%   c = criterion c (negative values -> bias towards yes responses)

% check input 
narginchk(2,2);

% check for values out of bounds, also issue errors or warnings if =1 or 0
if or(or(h>1,h<0),or(fA>1,fA<0))
    error('input arguments must fall in the 0 to 1 range')
elseif or(or(h==1,h==0),or(fA==1,fA==0))
    warning('This function will not return finite values when h or fA = 0 or 1.')
end

% d prime = z(h)-z(fA)
dp = norminv(h)-norminv(fA);

% c = -0.5*[z(h)+z(fA)]
c = -0.5*(norminv(h)+ norminv(fA));

end