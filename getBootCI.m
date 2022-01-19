function [ci] = getBootCI(vals, alpha)
%getBootCI Calculation of bootstrap confidence intervals for 'vals'
%times, 10000 bootstrap samples per default
%   'vals' are the input values. First dimension are observations
%   (e.g.trials): vals = [trials x time]

nBoots = 10;
bootType = 'per'; % 'per' or  'bca' (slow)

f  = @(x) mean(x); % bootstrap function is mean

rng(1,'twister'); % for reproducibility
ci = bootci(nBoots,{f,vals},'alpha',alpha,'type',bootType);

end
