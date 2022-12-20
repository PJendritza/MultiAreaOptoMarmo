function [pgh] = plotBootDistr(meanPl, CIs, xPlotPos, samples)
%plotBootDistr function to plot the distribution of all values from the boostrapping procedure
%   

% plotting parameters
alpha_pl = 0.001; % alpha for plotting range of distribution
distWidth = 0.6; % plotting width of the distribution

% create ks density from full bootstrapped data
[ksf, ksx] = ksdensity(samples);

% normalize ks density for plotting
ksf_norm = ksf/max(ksf)/2+xPlotPos;

% remove points outside the confidence intervals for plotting if required

CI_pl = prctile(samples,[100*alpha_pl*distWidth,100*(1-alpha_pl/2)]);

keepT = ksx>CI_pl(1) & ksx<CI_pl(2);
% keepT = ones(size(ksx));
ksf_norm_pl = ksf_norm(keepT);
ksx_pl = ksx(keepT);


hold on
pgon = polyshape(ksf_norm_pl, ksx_pl);
pgh = plot(pgon);
pgh.FaceColor = [0 0.6 1];
pgh.FaceAlpha = 0.4;
pgh.LineStyle = 'none';
plot([xPlotPos xPlotPos], [CIs(1) CIs(2)], 'k', 'LineWidth', 1)
plot(xPlotPos, meanPl, 'ko', 'MarkerSize',6, 'MarkerFaceColor','w')
hold off