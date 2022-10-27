function [sacc] = detectSaccEp(EYESeSc)
%detectSacc Detects saccades based on velocity in epoched data
%   Detailed explanation goes here
% EYESeSc: Offset corrected eye signal (x and y) in degress of visual angle

pxPerDeg = 27.854; % pixels per degree of visual angle at the center
Fs = 1000; % 1kHz sampling rate

monWidthD = 1680/pxPerDeg; % full width of the monitor
monHeighD = 1050/pxPerDeg; % full height of the monitor

% copy variable to process
tmp = EYESeSc;

% find points that are out of the monitor
outOfMonX = tmp(:,:,1)>monWidthD/2 | tmp(:,:,1)<-monWidthD/2;
outOfMonY = tmp(:,:,2)>monHeighD/2 | tmp(:,:,2)<-monHeighD/2;

% smooth the data with a gaussian window
msSmooth = 3;
tmp(:, :, 1) = gaussSmooth(tmp(:,:,1),Fs, msSmooth)';
tmp(:, :, 2) = gaussSmooth(tmp(:,:,2),Fs, msSmooth)';

% calculate eye velocity
EYEV(:, :, 1) = diff(tmp(:, :, 1),1, 2)*1000; % smooth velocity in deg/s
EYEV(:, :, 2) = diff(tmp(:, :, 2),1, 2)*1000; % smooth velocity in deg/s

EYEVabs = hypot(EYEV(:, :, 1), EYEV(:, :, 2)); % absulute  velocity
EYEAabs = diff(EYEVabs, 1, 2); % absulute  acceleration
EYEPabs = hypot(tmp(:, :, 1), tmp(:, :, 2)); % absulute  position

% blink detection based on jerk and out of monitor position
blinkJerkThresh = 100; % deg/s^3

isEYEblink = [ false(size(EYEVabs,1),3) abs(diff(EYEAabs(:,:), 1, 2))>blinkJerkThresh] | outOfMonX(:,:) | outOfMonY(:,:);

tCutBlinkPre = 25; % overwrite data -ms around blinks
tCutBlinkPost = 55; % overwrite data +ms around blinks
nTrls = size(EYESeSc,1);
nTpts = size(EYEVabs,2);

EYEVabs(isEYEblink(:,2:end)) = NaN;   % overwrite data +-25ms around blinks

for iTr = 1:nTrls
    % find blink on and off times
    blinkOnT = find(diff(isEYEblink(iTr, :))==1);
    blinkOffT = find(diff(isEYEblink(iTr, :))==-1);
    
    for iBon = 1:numel(blinkOnT)
        EYEVabs(iTr, max([blinkOnT(iBon)-tCutBlinkPre 1]):blinkOnT(iBon)) = NaN;
    end
    
    for iBoff = 1:numel(blinkOffT)
        EYEVabs(iTr, blinkOffT(iBoff):min([blinkOffT(iBoff)+tCutBlinkPost nTpts])) = NaN;
    end
end

% parameters for saccade detection
minPeakDist = 30; % ms
minPeakHeight = 50; % deg/s
maxPeakHeight = 1000; % deg/s
%%
disp('Finding saccades ... ')
for iTr = 1:nTrls
    [pksV{iTr}, lcsV{iTr}] = findpeaks(EYEVabs(iTr,:),'MinPeakDistance', round(minPeakDist*(Fs/1000)), 'MinPeakHeight', minPeakHeight);

    plotFlag = false; % manual flag to enable plotting for debug
    
    if plotFlag == true
        clf; hold on
        plot(tmp(iTr,:,1));
        plot(tmp(iTr,:,2));
        % plot(EYEPabs(iTr,:));
        axis tight
        plot(EYEVabs(iTr,:)/100, 'b');
        plot(lcsV{iTr}, pksV{iTr}/100, 'or')
    end
    
    medFbins = 50; % number of bins for the median filter
    
    % find saccade on and offset based on running median
    for i = 1:numel(pksV{iTr})
        
        trVData = EYEVabs(iTr,:);
        trVDataMed = medfilt1(EYEVabs(iTr,:),medFbins);
        
        thisSig = trVData(lcsV{iTr}(i):end)-trVDataMed(lcsV{iTr}(i):end);
        
        thisOffset = lcsV{iTr}(i)+find(thisSig <= 0 ,1,'first')+1; % add one sample for derivative difference
        if ~isempty(thisOffset) && pksV{iTr}(i) < maxPeakHeight
            saccOffsets{iTr}(i) = thisOffset;
        else
            saccOffsets{iTr}(i) = NaN; % in case onset was out of data
        end
        
        thisSig = trVData(1:lcsV{iTr}(i))-trVDataMed(1:lcsV{iTr}(i));
        
        thisOnset = lcsV{iTr}(i)-find(fliplr(thisSig) <= 0 ,1,'first')+1;% add one sample for derivative difference
        if ~isempty(thisOnset) && pksV{iTr}(i) < maxPeakHeight
            saccOnsets{iTr}(i) = thisOnset;
        else
            saccOnsets{iTr}(i) = NaN; % in case onset was out of data
        end
        
        % calculate saccade amplitude
        try
            xOn{iTr}(i)  = double(tmp(iTr,saccOnsets{iTr}(i),1));
            xOff{iTr}(i) = double(tmp(iTr,saccOffsets{iTr}(i),1));
            yOn{iTr}(i)  = double(tmp(iTr,saccOnsets{iTr}(i),2));
            yOff{iTr}(i) = double(tmp(iTr,saccOffsets{iTr}(i),2));
            
            saccAmpl{iTr}(i) = double(hypot(xOn{iTr}(i)-xOff{iTr}(i), yOn{iTr}(i)-yOff{iTr}(i)));
        catch
            xOn{iTr}(i)  = double(NaN);
            xOff{iTr}(i) = double(NaN);
            yOn{iTr}(i)  = double(NaN);
            yOff{iTr}(i) = double(NaN);
            
            saccAmpl{iTr}(i)   = double(NaN);
        end
        
        saccTrNr{iTr}(i) = iTr;
        if plotFlag == true
            try xline(saccOnsets{iTr}(i), 'k:'); end
            try xline(saccOffsets{iTr}(i), 'k:'); end
        end
    end
    
    if plotFlag == true
        axis tight
        ylim([-10 10])
        xlabel('Time (ms)')
        title('Saccade detecion')
    end
    
    if mod(iTr,25)==1
        disp(['Trial ' num2str(iTr) '/' num2str(nTrls)] )
    end
end

%% Create data vectors and structure 
sacc.pksV = pksV;
sacc.Ampl = saccAmpl;
sacc.Onsets = saccOnsets;
sacc.Offsets = saccOffsets;
sacc.TrNr = saccTrNr;
sacc.EYEVabs = EYEVabs;
sacc.xOn = xOn;
sacc.xOff = xOff;
sacc.yOn = yOn;
sacc.yOff = yOff;

allPeakVel = cell2mat(pksV);
allAmpl = cell2mat(saccAmpl);
allsaccOnsets = cell2mat(saccOnsets);
allsaccOffsets = cell2mat(saccOffsets);
allsaccTrNr = cell2mat(saccTrNr);

end

