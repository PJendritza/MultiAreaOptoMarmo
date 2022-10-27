function [ahOri] = oriTuningSubplot(tr, clu, iCluInd, trSel_yell, startBinT, endBinT)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here


% orientatin tuning analysis

spkCell = clu.spikeTimesEp(:,iCluInd);

dirConds = unique(tr.grat1direction);

% startBinT = 0.1;
% endBinT =  0.3;

binWidth = 0.001; % in seconds, bin to milliseconds

edges = [startBinT:binWidth:endBinT];
edgesCell = repmat({edges},size(spkCell,1),size(spkCell,2));

spks = cell2mat(cellfun(@histcounts,spkCell,edgesCell,'UniformOutput',false));

bootCIalpha = 0.05;

clear sua_cond ori ntr dur_ms sua_std sua_bootCI
for i = 1:numel(dirConds)
    
    thisDir = dirConds(i);
    [ntr(1,i), dur_ms(1,i)] =  size(spks(tr.grat1direction==thisDir & trSel_yell,:));
    sua_cond(1,i) =  sum(sum(spks(tr.grat1direction==thisDir & trSel_yell ,:)))/dur_ms(1,i)/ntr(1,i)*1000; % normalization to spks/s
    sua_std(1,i) =  std(sum(spks(tr.grat1direction==thisDir & trSel_yell ,:),2)./dur_ms(1,i).*1000);
    sua_sem(1,i) =  sem(sum(spks(tr.grat1direction==thisDir & trSel_yell ,:),2)./dur_ms(1,i).*1000);
    sua_bootCI(1:2,i) = getBootCI(sum(spks(tr.grat1direction==thisDir & trSel_yell ,:),2)./dur_ms(1,i).*1000, bootCIalpha);
    sua_cell{i} = sum(spks(tr.grat1direction==thisDir & trSel_yell ,:),2)./dur_ms(1,i).*1000;
    ori(1,i) = dirConds(i);
    
end

%%

%% von Mises fit

xaxM_rad = deg2rad(ori*2);
[thetaEst, kappaEst] = circ_vmpar(xaxM_rad',sua_cond');

opts_VM = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts_VM.Display = 'No';
% opts_Wedges.Startpoint = [0 180 0];
opts_VM.Lower = [0  -500 -pi];
opts_VM.Upper = [Inf 500 pi];
opts_VM.StartPoint=[1,10,thetaEst];


vonMisesFun = 'c*1/(2*pi*besseli(0,k))*exp(k*cos(x-m))'; % define von Mises function

normMua = max(sua_cond); %trapz(xaxWedgesVM_rad, aMUA_WedgesS); % get area to normalize data

[fitresultVM, gofVM] = fit( xaxM_rad', sua_cond'./normMua, vonMisesFun ,opts_VM);

fitPlotRangeVM = 0:0.1:2*pi;
FitOutVM = feval(fitresultVM,fitPlotRangeVM)*normMua;

subplot(4,3,12)
plot(rad2deg(fitPlotRangeVM/2),FitOutVM, 'k')
hold on

for i = 1:numel(ori)
    plot([ori(i) ori(i)],[sua_cond(i)-sua_sem(i) sua_cond(i)+sua_sem(i)], 'color', [0 0 0]);
    % plot([ori(i) ori(i)],[sua_cond(i)+sua_bootCI(2,i) sua_cond(i)-sua_bootCI(1,i)], 'color', [0 0 0])
    % scatter(ones(size(sua_cell{i}))*ori(i),sua_cell{i}, 10, 'filled','MarkerFaceColor',[0.5 0.5 0.5]', 'MarkerEdgeColor',[0 0 0],'jitter','on','jitterAmount',0.15);
end

plot(ori,sua_cond,'k.')

ylabel('spks/s')
% xlabel('Orientation (°)')
xlim([0 180])
box off
set(gcf,'Color','w')

ahOri = gca;
ahOri.XAxis.TickValues = [0 90 180];
xtickformat('%.0f°')
ahOri.Position = [0.5842    0.6650    0.0483    0.0463];
ahOri.FontSize = 6;

set(gca,'TickDir','out');

end

