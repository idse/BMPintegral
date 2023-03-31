clear all; close all;
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);

% raw data location, modify if not the same as the location of this script
dataDir = scriptPath;
dataDir = '/Users/idse/data_tmp/Emily/AW8 triplerep_ab_DAPI_stain_EF_190918';

meta = Metadata(dataDir);

%%
% manual metadata
%-------------------

meta.nWells = 4;
meta.posPerCondition = 5;

if meta.nPositions ~= meta.nWells*meta.posPerCondition
    warning('total position doesnt match nWells x posPerCondition');
end

% in order of imaging (clockwise, which is not the layout of the dish)
meta.conditions = {'C0','C1','C3','C10'};

% set this to true if making an '8-well' loop through a 4-well
loop4well = false;

meta.nucChannel = 1;

% name your channels here if you want
BraChannel = 0;
nucChannel = 1;

%% save stitched previews of the MIPs

stitchedPreviews(dataDir, meta); 

%% extract nuclear and cytoplasmic levels

% notes:
% Andor channel counting starts at 0
% with cytoplasm: 'cytoplasmicLevels',true, 'cytoSize', 8 (example)

opts = struct(...
                    'dataChannels',     [BraChannel nucChannel],...
                    'fgChannel',        BraChannel,...
                    'cytoplasmicLevels',false,... 
                    'segmentationDir',  fullfile(dataDir,'MIP'),...
                    'nucShrinkage',     2,...
                    'bgMargin',         10,...
                    'segFG',            0);

opts.cleanupOptions = struct('separateFused', true,...
    'clearBorder',true, 'minAreaStd', 1, 'minSolidity',0, 'minArea', 150);


%% check that the options are set right

pi = 1;
P = Position(meta);
P.setID(pi);
ti = 1;
opts.tMax = ti;

dbInfo = P.extractData(dataDir, opts);

bg = P.cellData(ti).background
nucl = P.cellData(ti).nucLevelAvg
% cytl = P.cellData(ti).cytLevelAvg
% (nucl-bg)./(cytl - bg)

im = P.loadImage(dataDir, nucChannel, ti);
MIP = max(im,[],3);
A = 0.8*imadjust(mat2gray(MIP));
s = 0.3;
imshow(cat(3, A + dbInfo.nucmaskraw*s, A + s*dbInfo.nucmask, A + s*dbInfo.cytmask));

% try out the nuclear cleanup settings on some frame:
bla = nuclearCleanup(seg(:,:,time), opts.cleanupOptions);
imshow(bla)

%%

P.makeTracks(20);

figure,
s = 0.4;
% THERE MAY BE AN ISSUE that objectid=1 and misdetections are confused
% right from the Ilastik export
imshow(cat(3, A + dbInfo.nucmask, A + s*dbInfo.nucmaskraw, A + s*cytmask));

opts.tMax = meta.nTime;

% plot trajectories on image simpletracker

tmax = 100;
colors = jet(tmax);
XY = P.timeTraces.trackXY;
for pi = 1:numel(P.timeTraces.trackXY)

    x = squeeze(XY{pi}(:,1));
    y = squeeze(XY{pi}(:,2));
    p = line(x, y, 'LineWidth',2);
    colorsP = colors(~isnan(x),:);
    colorsP = [255*colorsP colorsP(:,1)*0+255];
    colorsP = uint8(colorsP');
    drawnow
    set(p.Edge, 'ColorBinding','interpolated', 'ColorData',colorsP)
end
hold off

%% plot trajectories on image

figure,
v = VideoWriter(fullfile(dataDir,'trackTest2.avi'));%,'Uncompressed AVI');
v.FrameRate = 5;
open(v)

tmax = 100;
for ti = 1:tmax
    
    im = P.loadImage(dataDir, BraChannel, ti);
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

    colors = jet(tmax);
    XY = P.timeTraces.trackXY;
    T = P.timeTraces.trackT;
    
    for pi = 1:numel(P.timeTraces.trackXY)

        idx = T{pi} <= ti;
        cidx = 0*(1:tmax);
        cidx(T{pi}) = T{pi};
        cidx(cidx > ti) = 0;
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

%% plot trajectories on image without running extractData

info = h5info(fullfile(MIPdir,'Tracking.h5'));
trackingLabels = h5read(fullfile(MIPdir,'Tracking.h5'),'/exported_data');
trackingLabels = squeeze(trackingLabels);

%%
CMall = [];
for ti = 1%:size(trackingLabels,3)
   
    trackL = trackingLabels(:,:,ti)';
    L = bwlabel(trackingLabels(:,:,ti)' > 1);
    props = regionprops(trackL, trackL, 'Centroid', 'MeanIntensity');
    CM = cat(1,props.Centroid);
    lineage = round([props.MeanIntensity]);
    CMall = cat(3,CMall,CM);
end

%%
figure,
v = VideoWriter(fullfile(dataDir,'trackTestX.avi'));%,'Uncompressed AVI');
v.FrameRate = 5;
open(v)

tmax = 100;
for ti = 1:tmax
    
    imshow(trackingLabels(:,:,ti)',[]);
    hold on

    colors = jet(tmax);
    XY = CMall;
    for pi = 2:46

        x = squeeze(XY(pi,2,1:ti));
        y = squeeze(XY(pi,1,1:ti));
        p = line(x, y, 'LineWidth',2);
        colorsP = colors(~isnan(x),:);
        colorsP = [255*colorsP colorsP(:,1)*0+255];
        colorsP = uint8(colorsP');
        drawnow
        set(p.Edge, 'ColorBinding','interpolated', 'ColorData',colorsP)
    end
    hold off

    writeVideo(v,getframe)
end
close(v);


%% all signaling traces

figure,
channel = 1;

bg = cat(3,P.cellData.background);
bg = squeeze(bg(:,channel,:));
N = cat(3,P.cellData.nucLevel);
N = squeeze(N(:,channel,:))';

C = cat(3,P.cellData.cytLevel);
C = squeeze(C(:,channel,:))';

plot((N-bg)./(C-bg))

%% some specific traces

positions = [10 20 30 40];
N = cat(3,P.cellData.nucLevel);
C = cat(3,P.cellData.cytLevel);
clf 
hold on
for pos = positions
    n = squeeze(N(pos,channel,:));
    c = squeeze(C(pos,channel,:));
    plot(1:100, (n-bg)./(c-bg));
end
hold off

%%
fname = 'Gridse_MIP_p0000_w0001_Object-Identities.h5';
fname = 'Gridse_MIP_p0000_w0001_Tracking-Result.h5';
info = h5info(fullfile(MIPdir,fname));
bla = h5read(fullfile(MIPdir,fname),'/exported_data');
ti = 5;
figure,imshow(squeeze(bla(1,:,:,ti))',[])

%%
fname = 'Gridse_MIP_p0000_w0001_H5-Event-Sequence.h5/00005.h5';
info = h5info(fullfile(MIPdir,fname));
blaa = h5read(fullfile(MIPdir,fname), '/segmentation/labels');
imshow(squeeze(blaa),[])



%% run the analysis on all time points

tic
positions(meta.nPositions) = DynamicPositionAndor();

for pi = 1%meta.nPositions

    positions(pi) = DynamicPositionAndor(meta, pi);
    positions(pi).extractData(dataDir, nucChannel, opts);
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

%% make figure of the time traces

options = struct('wellnrs', 1, 'channel', BraChannel,...
                'mode','N:C',...
                'yLabel', 'Smad4 (N:C)',...
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

%%
% multiple conditions 
options.wellnrs = wellnrs;
makeSignalingGraphs(positions, meta, options)
legend(meta.conditions);
export_fig(fullfile(dataDir,'timeTrace_multipleConditions.png'),'-native -m2');

