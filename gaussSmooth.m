function [ dataSmooth ] = gaussSmooth( data, Fs, sigma_ms )
%GAUSSSMOOTH Smoothes inputdata with sampling rate Fs along first dimention with a gaussian
%kernel of sigma milliseconds
%  
% Mar 2021: modified to use smoothdata.m, P. Jendritza

sigma_s = sigma_ms/1000; % one sigma in seconds
gCut = 2.5; % gauss window is cut at 2.5 sigma
gWinSize = round(Fs*sigma_s*gCut*2);

% w = gausswin(gWinSize, gCut);
% w = w/sum(w); % normalize window

dataSmooth = double(data)'; % allocate
% dataSmooth = filtfilt(w, 1, double(data)');

dataSmooth = smoothdata(data','gaussian',gWinSize);

end

