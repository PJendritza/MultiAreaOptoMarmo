% Figure 5_be
% Code to generate receptive field plots based on MUA from area V1 and V6

% load data
addpath data
load('Fig5_be_data.mat')

% % get channel selection for modulated channels
selV1chns = chInfo.isModuTBothV1;
selV6chns = chInfo.isModuTBothV6;


%% create data for plotting

figure('Name',['Fig5be'], 'pos',[10   10   945   454])

% % plot RFs for the selected shank
% addpath /mnt/hpx/home/jendritzap/MarmosetV1V6analysis/code/
% addpath /mnt/hpx/home/jendritzap/MarmosetV1V6analysis/code/ploting

imThresh = 0.2; % threshold for RF outline

% pmap = flipud(getMappingV1V6('V1_1A','normal', cfg.sessionID));
pmapV1 = [45;47;43;33;9;3;37;35;5;7;39;41;11;13;1;15;27;25;57;55;21;19;31;17;51;49;53;63;23;29;59;61];% channel mapping in V1

nREspCh = sum(selV1chns(pmapV1)); % number of resp. channels on this shank
outCol = gray(round(nREspCh*1.5)); % gray colors for ouline RF
outCol = flipud(outCol(1:nREspCh,:));

clear S
for i = 1:32
    ich = pmapV1(i);
    S(:,:,i) =    dataV1Cell{ich}.RFimgModel1;
    [RFcenterY(i), RFcenterX(i)] = pol2cart(-dataV1Cell{ich}.RF_polangVM,dataV1Cell{ich}.RF_eccent);
end


% axis scaling
xax = aux.xaxPx/aux.pxPerDeg;
yax = aux.yaxPx/aux.pxPerDeg;

[xaxStack,yaxStack] = meshgrid(xax,yax);
% Define the height Z of the surface at the coordinates given by (X,Y).
Zstack = ones(size(S(:,:,1)));
% Warp the image over the surface defined by the coordinates (X,Y,Z).

clf
subplot(1,2,1)


%% V1 plot
ip = 1;
for i = 1:32
    
    if selV1chns(pmapV1(i))
        
        
        iSimg = rescale(S(:,:,i));
        [~, maxInd] = max(iSimg(:));
        [peakY(i),peakX(i)]= ind2sub(size(iSimg),maxInd);
        
        h{i} = surf(xaxStack,yaxStack,Zstack*i, iSimg,'EdgeColor','none');
        hold on
                
        alphaMapSurf{i} = iSimg>imThresh;
        alphaMapSurf{i} = imgaussfilt(double(alphaMapSurf{i}), 1);
        
        set(h{i}, 'AlphaData', alphaMapSurf{i}, 'AlphaDataMapping', 'none');
        set(h{i}, 'FaceAlpha', 'flat');
        
        % plot RF outline
        iSimgBW = iSimg>imThresh;
        [B,L] = bwboundaries(iSimgBW,'noholes');
        iSimgBoundary = B{1};
        p = 1:numel(iSimgBoundary(:,1));
        clear iSimgBoundaryInterp
        xRFoutl = xax(iSimgBoundary(:,2));
        yRFoutl = yax(iSimgBoundary(:,1));
        iSimgBoundaryInterp(:,2) = smooth(interp1(p,xRFoutl,1:0.1:numel(iSimgBoundary(:,1)),'spline'),30); % inperpolate and smooth to reduce the pixel artifacts
        iSimgBoundaryInterp(:,1) = smooth(interp1(p,yRFoutl,1:0.1:numel(iSimgBoundary(:,2)),'spline'),30); % inperpolate and smooth to reduce the pixel artifacts
        plot(+iSimgBoundaryInterp(:,2), +iSimgBoundaryInterp(:,1), 'color', outCol(ip,:), 'LineWidth', 1)
        ip = ip+1;
    end
    
end

centerX = xax(median(peakX));
centerY = yax(median(peakY));

plot3([centerX centerX], [centerY centerY], [0 33],'k','linewidth',1)
plot3([-1 1], [0 0], [0 0],'k','linewidth',1)
plot3([0 0], [-1 1], [0 0],'k','linewidth',1)
xysk = 5;
xlim([-1 1]*xysk+centerX)
ylim([-1 1]*xysk+centerY)
zlim([0.0 32.5])
pbaspect([1 1 4])
zticks = 1:32;
zticks = zticks(1:2:end);
set(gca,'ztick', zticks)

ax= gca; fontSi = 8;
ax.XAxis.FontSize = fontSi;
ax.YAxis.FontSize = fontSi;
ax.ZAxis.FontSize = fontSi;

colormap(inferno)

xlabel('X-pos. (°)')
ylabel('Y-pos. (°)')
zlabel('Electrode site')
set(gcf,'color','w');

drawnow
disp('V1 plot done')


%% V6 plot
subplot(1,2,2)
cla

% pmap = fliplr(getMappingV1V6('V6_1A','normal', cfg.sessionID));
pmapV6 = [61;59;29;23;63;53;49;51;17;31;19;21;55;57;25;27;15;1;13;11;41;39;7;5;35;37;3;9;33;43;47;45];% channel mapping in V6

nREspCh = sum(selV6chns(pmapV6)); % number of resp. channels on this shank
outCol = gray(round(nREspCh*1.5)); % gray colors for ouline RF
outCol = flipud(outCol(1:nREspCh,:));

lowerZ = 5.5; % lower boundary for plot

clear S
for i = 1:32
    ich = pmapV6(i);
    S(:,:,i) =    dataV6Cell{ich}.RFimgModel6;
    
end

peakY = nan(1, 32);
peakX = peakY;

ip = 1;
for i = 1:32
    
    if selV6chns(pmapV6(i))
        iSimg = rescale(S(:,:,i));
        
        [~, maxInd] = max(iSimg(:));
        [peakY(i),peakX(i)]= ind2sub(size(iSimg),maxInd);
        
        h{i} = surf(xaxStack,yaxStack,Zstack*i, iSimg,'EdgeColor','none');
        hold on
        
        alphaMapSurf{i} = iSimg>imThresh;
        alphaMapSurf{i} = imgaussfilt(double(alphaMapSurf{i}), 3);
        
        set(h{i}, 'AlphaData', alphaMapSurf{i}, 'AlphaDataMapping', 'none');
        set(h{i}, 'FaceAlpha', 'flat');
        
        % plot RF outline
        iSimgBW = iSimg>imThresh;
        [B,L] = bwboundaries(iSimgBW,'noholes');
        iSimgBoundary = B{1};
        p = 1:numel(iSimgBoundary(:,1));
        clear iSimgBoundaryInterp
        xRFoutl = xax(iSimgBoundary(:,2));
        yRFoutl = yax(iSimgBoundary(:,1));
        iSimgBoundaryInterp(:,2) = smooth(interp1(p,xRFoutl,1:0.1:numel(iSimgBoundary(:,1)),'spline'),30); % inperpolate and smooth to reduce the pixel artifacts
        iSimgBoundaryInterp(:,1) = smooth(interp1(p,yRFoutl,1:0.1:numel(iSimgBoundary(:,2)),'spline'),30); % inperpolate and smooth to reduce the pixel artifacts
        plot3(+iSimgBoundaryInterp(:,2), +iSimgBoundaryInterp(:,1), ones(size(iSimgBoundaryInterp(:,1)))*lowerZ, 'color', outCol(ip,:), 'LineWidth', 1)
        ip = ip+1;
    end
    
end

centerX = xax(nanmedian(peakX));
centerY = yax(nanmedian(peakY));

plot3([centerX centerX], [centerY centerY], [lowerZ 33],'k','linewidth',1)
plot3([-4 4], [0 0], [1 1]*lowerZ,'k','linewidth',1)
plot3([0 0], [-4 4], [1 1]*lowerZ,'k','linewidth',1)

xysk = 23;
xlim([-1 1]*xysk+centerX)
ylim([-1 1]*xysk+centerY)

zlim([lowerZ 32.5])
pbaspect([1 1 4])
zticks = 1:32;
zticks = zticks(1:2:end);
set(gca,'ztick', zticks)

ax= gca; fontSi = 8;
ax.XAxis.FontSize = fontSi;
ax.YAxis.FontSize = fontSi;
ax.ZAxis.FontSize = fontSi;

colormap(inferno)

xlabel('X-pos. (°)')
ylabel('Y-pos. (°)')
zlabel('Electrode site')
set(gcf,'color','w');

drawnow
disp('V6 plot done')

