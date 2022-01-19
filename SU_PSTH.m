function [ph] = SU_PSTH(spkCell,msSmooth, startBinT, endBinT)
%plotSpikePSTH Function to plot smooth PSTH from spike time data
%   Detailed explanation goes here

binWidth = 0.001; % in seconds, bin to milliseconds
edges = [startBinT:binWidth:endBinT];

%
% tmp = cell2mat(spkCell);
% nTr = size(spkCell,1);
% [N,~] = histcounts(tmp,edges);

edgesCell = repmat({edges},size(spkCell,1),size(spkCell,2));

spks = cell2mat(cellfun(@histcounts,spkCell,edgesCell,'UniformOutput',false));
spks = spks/binWidth; % normalize to seconds with binwidth

%% smoothing with Gaussian window
Fs = 1/binWidth;

spks_smooth = gaussSmooth( spks, Fs, msSmooth)';

%% calculate SEM
if size(spks_smooth,1) ==1
    plotErr = false;
else
    spks_smoothSEM = sem(spks_smooth);
    plotErr = true;
end

%% create time vector
edgewidth = diff(edges(1:2));
timeX = edges(1:end-1)+edgewidth/2;

%% plot shaded error bar

if plotErr
    ph = shadedErrorBar(timeX,mean(spks_smooth), spks_smoothSEM,'lineprops',{'color', 'k','linewidth',0.5},'transparent',1);
else
    ph = plot(timeX,spks_smooth, 'color', 'k','linewidth',0.5);
end

