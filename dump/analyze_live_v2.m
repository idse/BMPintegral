clear all; close all;
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);

% raw data location, modify if not the same as the location of this script
dataDir = scriptPath;
MIPdir = fullfile(dataDir,'MIP');
cd(dataDir);

meta = MetadataAndor(dataDir);

% manual metadata
%-------------------

treatmentTime = 4;

meta.nWells = 1;
meta.posPerCondition = 1;

if meta.nPositions ~= meta.nWells*meta.posPerCondition
    warning('total position doesnt match nWells x posPerCondition');
end

% in order of imaging (clockwise, which is not the layout of the dish)
meta.conditions = {'?'};

% set this to true if making an '8-well' loop through a 4-well
loop4well = false;

nucChannel = 1;
meta.nucChannel = 1;
S4Channel = 0;

tres = 5/60;

%% save stitched previews of the MIPs

stitchedPreviews(dataDir, meta); 

%% extract nuclear and cytoplasmic levels

% Andor channel counting starts at 0
opts = struct(...
                    'dataChannels',     [S4Channel nucChannel],...
                    'fgChannel',        S4Channel,...
                    'cytoplasmicLevels',true,... 
                    'segmentationDir',  fullfile(dataDir,'MIP'),...
                    'nucShrinkage',     2,...
                    'cytoSize',         8,...
                    'bgMargin',         10,...
                    'NCRcutoff',        [3 Inf],...
                    'segFG',            0);

opts.cleanupOptions = struct('separateFused', true,...
    'clearBorder',true, 'minAreaStd', 1, 'minSolidity',0, 'minArea',500);


%% check that the options are set right

pi = 1;
P = DynamicPositionAndor(meta, pi);
P.setID(pi);
ti = 100;
opts.tMax = ti;

% try out the nuclear cleanup settings on some frame:
% bla = nuclearCleanup(seg(:,:,time), opts.cleanupOptions);
% imshow(bla)
% 
% dbInfo = P.extractData(dataDir, opts);
% 
% cytmask = false(size(dbInfo.nucmask));
% cytmask(cat(1,dbInfo.cytCC.PixelIdxList{:}))=true;

dbInfo = P.extractData(dataDir, opts);
cytmask = dbInfo.Lcyt > 0;

bg = P.cellData(ti).background
nucl = P.cellData(ti).nucLevelAvg
cytl = P.cellData(ti).cytLevelAvg
(nucl-bg)./(cytl - bg)

im = P.loadImage(dataDir, S4Channel, ti);
%labels = P.loadTrackingLabels(MIPdir, nucChannel, time);
MIP = max(im,[],3);
A = imadjust(mat2gray(MIP));

%% run the analysis on all time points

tic
positions(meta.nPositions) = DynamicPositionAndor();

for pi = 1%meta.nPositions

    positions(pi) = DynamicPositionAndor(meta, pi);
    positions(pi).extractData(dataDir, opts);
    positions(pi).makeTimeTraces();
    save(fullfile(dataDir,'positions'), 'positions');
end
toc

%% load results if above block was run previously

load(fullfile(dataDir,'positions'));

% % to concatenate two files (if microscope had to be stopped)
%
% positions_before = load(fullfile(dataDir,'positions_before'));
% positions_before = positions_before.positions;
% 
% for pi = 1:meta.nPositions
%     positions(pi).cellData = cat(2, positions_before(pi).cellData,...
%                                             positions(pi).cellData);
%     positions(pi).ncells = cat(2, positions_before(pi).ncells,...
%                                             positions(pi).ncells);
% 	positions(pi).nTime = positions_before(pi).nTime + positions(pi).nTime;
%     positions(pi).makeTimeTraces;
% end
% 
% save(fullfile(dataDir,'positions'), 'positions');

%% make tracks

P = positions(1);
P.makeTracks(20);

nTracks = numel(P.timeTraces.trackXY);

trackMatrix = NaN([meta.nTime nTracks 2]);

for tri=1:nTracks 
    XY = P.timeTraces.trackXY;
    T = P.timeTraces.trackT;
    trackMatrix(T{tri}, tri, :) = XY{tri};
end

disp(['longest track: ' num2str(max(cellfun(@numel,P.timeTraces.trackT)))]);

% load segmentation
seg = P.loadSegmentation(MIPdir, nucChannel);


%% plot signaling traces

nTracks= numel(P.timeTraces.nucLevel);
N = NaN([tmax nTracks]);
C = NaN([tmax nTracks]);
B = cat(1,P.cellData.background);

ci = 1;
for tri = 1:nTracks
    N(P.timeTraces.trackT{tri}, tri) = P.timeTraces.nucLevel{tri}(:,ci);
    C(P.timeTraces.trackT{tri}, tri) = P.timeTraces.cytLevel{tri}(:,ci);
end
B = B(:,ci);
traceMatrix = (N-B)./(C-B);

% all signaling traces
time = (1:tmax)*tres;
plot(time, traceMatrix)
ylim([0.5 2]);
xlim([0 70]);
set(gca, 'LineWidth', 2);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')
xlabel('time (hrs)');
ylabel('signaling');
saveas(gcf,'allTimeTraces.png');

% % some specific traces
% figure,
% mytracks = find(longestidx);
% S = (N(:,mytracks)-B)./(C(:,mytracks)-B);
% 
% hold on
% for i = 1:numel(mytracks)
%     
%     ti = find(~isnan(S(:,i)));
%     time = ti*5/60;
%     plot(time, S(~isnan(S(:,i)),i));
% end
% hold off


%% PCA on signaling traces

traceMatrix(isnan(traceMatrix))=0;
S = cov(traceMatrix');
[V,D] = eigs(S);
plot(V(:,1)')
hold on 
plot(V(:,2)')
plot(V(:,3)')
hold off


%% plot some tracks on the final frame

figure,

ti = meta.nTime;

im = P.loadImage(dataDir, S4Channel, ti);
MIP = max(im,[],3);
A = imadjust(mat2gray(MIP));

nucmask = seg(:,:,ti);
nucmask = nucmask - imerode(nucmask,strel('disk',3));
imshow(cat(3, A + 0*nucmask, A + nucmask, A));

longestidx = cellfun(@numel,T) > 600;
hold on 

tri = find(longestidx);
trackColors = lines(numel(tri));

for i = 1:numel(tri)
    X = trackMatrix(:,tri(i),1);
    Y = trackMatrix(:,tri(i),2);
    p = line(X(~isnan(X)), Y(~isnan(X)), 'LineWidth',2,'Color',trackColors(i,:));
end

%p = line(trackMatrix(:,longestidx,1), trackMatrix(:,longestidx,2), 'LineWidth',2);
%scatter(trackMatrix(end,longestidx,1), trackMatrix(end,longestidx,2),'filled')
hold off

% for pi = 1:numel(P.timeTraces.trackXY)
% 
%     x = squeeze(XY{pi}(:,1));
%     y = squeeze(XY{pi}(:,2));
%     p = line(x, y, 'LineWidth',2);
%     colorsP = colors(~isnan(x),:);
%     colorsP = [255*colorsP colorsP(:,1)*0+255];
%     colorsP = uint8(colorsP');
%     drawnow
%     set(p.Edge, 'ColorBinding','interpolated', 'ColorData',colorsP)
% end
% hold off
% 

%% plot trajectories on video

% P = positions(1);
% P.makeTracks(20);

v = VideoWriter(fullfile(dataDir,'trackTest5.avi'));%,'Uncompressed AVI');
v.FrameRate = 5;
open(v)

trackColors = jet(nTracks);
trackColors = trackColors(randperm(nTracks),:);
    

figure,
tmax = meta.nTime;
for ti = 2:tmax
    
    fprintf('.');
    if mod(ti,100)==0
        fprintf('\n');
    end
    im = P.loadImage(dataDir, S4Channel, ti);
    MIP = max(im,[],3);
    A = mat2gray(MIP);
    A = imadjust(A, stretchlim(A, 0.001));
    nucmask = seg(:,:,ti);
    nucmask = nucmask - imerode(nucmask,strel('disk',3));
    imshow(cat(3, A, A + nucmask, A));

    hold on

    %scatter(trackMatrix(ti,:,1),trackMatrix(ti,:,2),[],trackColors,'filled');
    
    longestidx = cellfun(@numel,T) > 600;

    %tri = 1:nTracks;
    tri = longestidx;
    p = line(trackMatrix(1:ti,tri,1), trackMatrix(1:ti,tri,2), 'LineWidth',2);
    
    hold off

    writeVideo(v,getframe)
end
close(v);

%% video of trajectories and signaling traces

v = VideoWriter(fullfile(dataDir,'trackTest7.avi'));%,'Uncompressed AVI');
v.FrameRate = 15;
open(v)

margins = [0.06,0.02];
n = 4; m = 7;

longestidx = cellfun(@numel,T) > 600;
%tri = 1:nTracks;
tri = find(longestidx);
S = (N(:,tri)-B)./(C(:,tri)-B);

%trackColors = jet(nTracks);
%trackColors = trackColors(randperm(nTracks),:);
tri = tri(1:6);
trackColors = lines(6);

F = figure('Position',[1 1 2000 1000]);

for ti = 2:tmax
    
    clf
    
    fprintf('.');
    if mod(ti,100)==0
        fprintf('\n');
    end
    
    im = P.loadImage(dataDir, S4Channel, ti);
    MIP = max(im,[],3);
    A = mat2gray(MIP);
    A = imadjust(A, stretchlim(A, 0.001));
    nucmask = seg(:,:,ti);
    nucmask = nucmask - imerode(nucmask,strel('disk',3));
    
    subplot_tight(n,m,[1:3 m+(1:3) 2*m+(1:3)],margins);
    imshow(cat(3, A, A + nucmask, A));

    hold on

    %scatter(trackMatrix(ti,:,1),trackMatrix(ti,:,2),[],trackColors,'filled');
    
    for i = 1:numel(tri)
        X = trackMatrix(1:ti,tri(i),1);
        Y = trackMatrix(1:ti,tri(i),2);
        p = line(X(~isnan(X)), Y(~isnan(X)), 'LineWidth',2,'Color',trackColors(i,:));
    end
    %h = plot(trackMatrix(1:ti,tri,1), trackMatrix(1:ti,tri,2), 'LineWidth',2);
    %set(h, {'color'}, num2cell(trackColors,2));%trackColors(tri,:)
    %set(p.Edge, 'ColorBinding','interpolated', 'ColorData',colorsP)
    
    fs = 24;
    set(gcf,'color','w');
    
    for I = 1:3
    
        subplot_tight(n,m, m*I-1:m*I,margins);
        timeidx = find(~isnan(S(1:ti,1)));
        time = timeidx*tres;
        plot(time, S(timeidx, I),'Color',trackColors(I,:),'LineWidth',2);
        axis([0 meta.nTime*tres 0.5 1.5]);
        set(gca, 'LineWidth', 2);
        set(gca,'FontSize', fs)
        set(gca,'FontWeight', 'bold')
        yticklabels([]);
    end
    xlabel('time (hrs)');
    for I = 4:6
    
        subplot_tight(n,m, m*(I-3)-3:m*(I-3)-2,margins);
        timeidx = find(~isnan(S(1:ti,1)));
        time = timeidx*tres;
        plot(time, S(timeidx, I),'Color',trackColors(I,:),'LineWidth',2);
        axis([0 meta.nTime*tres 0.5 1.5]);
        set(gca, 'LineWidth', 2);
        set(gca,'FontSize', fs)
        set(gca,'FontWeight', 'bold')
        
        if I == 5
            ylabel('signaling');
        end
    end
    xlabel('time (hrs)');

    hold off

    writeVideo(v,getframe(F))
end

close(v);

%% plot trajectories on video with temporal color gradient

figure,
v = VideoWriter(fullfile(dataDir,'trackTest3.avi'));%,'Uncompressed AVI');
v.FrameRate = 5;
open(v)

tmax = meta.nTime;
for ti = 1:tmax
    
    fprintf('.');
    if mod(ti,100)==0
        fprintf('\n');
    end
    im = P.loadImage(dataDir, S4Channel, ti);
    %labels = P.loadTrackingLabels(MIPdir, nucChannel, time);
    MIP = max(im,[],3);
    A = imadjust(mat2gray(MIP));
    nucmask = s*dbInfo.segall{2}(:,:,ti);
    nucmask = nucmask - imerode(nucmask,strel('disk',3));
    imshow(cat(3, A + 0*dbInfo.nucmask, A + nucmask, A));
    %trackLabels = this.loadTrackingResult(opts.segmentationDir);
    %imshow(trackLabels(:,:,ti),[]);
    %imshow(label2rgb(trackLabels(:,:,100),'gray','k','shuffle'));
    hold on

    twindow = 100;
    colors = jet(tmax);
    XY = P.timeTraces.trackXY;
    T = P.timeTraces.trackT;
    
    for pi = 1:numel(P.timeTraces.trackXY)

        idx = T{pi} <= ti & T{pi} > ti - twindow;
        cidx = 0*(1:tmax);
        cidx(T{pi}) = T{pi};
        cidx(cidx > ti) = 0;
        cidx(cidx <= ti - twindow) = 0;
        cidx = cidx > 0;
        x = squeeze(XY{pi}(:,1));
        y = squeeze(XY{pi}(:,2));
        if sum(idx)>0
            p = line(x(idx), y(idx), 'LineWidth',2);
            colorsP = colors(cidx,:);
            colorsP = [255*colorsP colorsP(:,1)*0+255];
            colorsP = uint8(colorsP');
            drawnow
            set(p.Edge, 'ColorBinding','interpolated', 'ColorData',colorsP)
        end
    end
    hold off

    writeVideo(v,getframe)
end
close(v);


%% make figure of the time traces

options = struct('wellnrs', 1, 'channel', S4Channel,...
                'loop4well', loop4well,...
                'minNCells', 10,...
                'dataChannels', opts.dataChannels,...
                'treatmentTime', treatmentTime,...
                'signalingRange', [0.6 1.2]);

% individual conditions
wellnrs = options.wellnrs;
for wellnr = wellnrs
    options.wellnrs = wellnr;
    makeSignalingGraphs(positions, meta, options);
    export_fig(fullfile(dataDir,['timeTrace_well' num2str(wellnr) '.png']),'-native -m2');
end

% multiple conditions 
options.wellnrs = wellnrs;
makeSignalingGraphs(positions, meta, options)
legend(meta.conditions);
export_fig(fullfile(dataDir,'timeTrace_multipleConditions.png'),'-native -m2');

