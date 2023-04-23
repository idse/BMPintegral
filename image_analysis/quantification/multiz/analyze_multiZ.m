clear; close all; clc

%% Set options
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
dataDir = scriptPath;
imdir = dataDir;
bare = 'stitched_p%.4d_w%.4d_t%.4d';
fnameformat = [bare,'.tif'];

manualMeta = struct();

% number of experimental conditions
manualMeta.nWells = 4;
manualMeta.posPerCondition = 5;
manualMeta.nPositions = manualMeta.nWells*manualMeta.posPerCondition;
manualMeta.plateLayoutNames = false;
% set this to true if making an '8-well' loop through a 4-well
manualMeta.loop4well = false;
% as plate layout
manualMeta.conditions = {'160k,BMP4+IWP2','160k,BMP4',...
    '240k,BMP4','240k,BMP4+IWP2'};

manualMeta.nucChannel = 0;
manualMeta.channelLabel = {'DAPI','ISL1','pSmad1','TBXT'};
manualMeta.timeInterval = '20 min';
% complete the rest of the metadata automatically and make metadata object
meta = Metadata(dataDir, manualMeta);
save(fullfile(dataDir,'meta'), 'meta');
% name your channels here if you want
nucChannel = meta.nucChannel;
smadChannel = 1;
% specify if timelapse:
treatmentTime = 3; % frame of treatment

timeres = strsplit(meta.timeInterval,' ');
timeres = str2double(timeres{1});
tscale = timeres/60;
ntime = meta.nTime;
tvec = (-treatmentTime:ntime - 1 - treatmentTime)'*tscale;

%% cleanup and quantification options
opts = struct(...
                    'dataChannels',     0:3,...
                    'nucChannel',       nucChannel,...
                    'cytoplasmicLevels',true,...
                    'fgChannel',        [],...
                    'cytoSize',         4,...
                    'cytoMargin',       0,...
                    'segmentationDir',  dataDir,...
                    'nucShrinkage',     2,...
                    'bgMargin',         4);

zopts = struct('IoU',0.6,'maxZsize',15,'maxCentroidDistance',15,'Zmethod','useNuclearMarker');
segdir = opts.segmentationDir;

%% Read in a single image and quantify the nuclei + cytoplasm in it
close all
pidx = 11;
P = Colony(meta,pidx);
P.setRadius(350, meta.xres);
P.dataChannels = 1:P.nChannels;
P.filenameFormat = fnameformat;
ti = 1;

seg = cell(1,meta.nChannels);
img = cell(1,meta.nChannels);
for cii = 1:meta.nChannels
    ci = opts.dataChannels(cii);
    img{ci+1} = P.loadImage(imdir,ci,ti);
end

%read nuclear segmentation
name = fullfile(segdir,sprintf(bare,pidx-1,nucChannel,ti-1));
segname = [name,'_FinalSegmentation.tif'];
nucseg = false(size(img{1},1),size(img{1},2),meta.nZslices);
for zi = 1:meta.nZslices
    nucseg(:,:,zi) = imread(segname,zi);
end
seg{nucChannel+1} = nucseg;

%read cytoplasmic segmentation
if ~isempty(opts.fgChannel)
    name = fullfile(segdir,sprintf(bare,pidx-1,opts.fgChannel,ti-1));
    segname = [name,'_Simple Segmentation.h5'];
    fgseg = ilastikRead(segname);
    seg{opts.fgChannel+1} = fgseg;
end

cellData = struct;
figure('WindowState','maximized')
for zi = 1:meta.nZslices
    ims = cellfun(@(x) x(:,:,zi),img,'UniformOutput',false);
    segs = cell(1,meta.nChannels);
    segs{nucChannel+1} = seg{nucChannel+1}(:,:,zi);
    if ~isempty(opts.fgChannel)
        segs{opts.fgChannel+1} = seg{opts.fgChannel+1}(:,:,zi);
    end
    
    opts.nuclearSegmentation = segs{nucChannel+1};
    [dbInfo, stats] = extractDataSingleFrame(ims,segs,opts);
    fields = fieldnames(stats);
    for fi = 1:length(fields)
        cellData(zi).(fields{fi}) = stats.(fields{fi});
    end
end

[newCellData, chain] = combineZsegmentation(cellData, meta, zopts);
% [newCellData, chain] = combineZtest2(cellData, meta, zopts);
fields = fieldnames(newCellData);
for fi = 1:length(fields)
    P.cellData(ti).(fields{fi}) = newCellData.(fields{fi});
end
P.ncells(ti) = size(P.cellData(ti).XY,1);

% P.setCenter();
% P.makeRadialAvgSeg();

%troubleshooting
XY = newCellData.XY;
ax = gobjects(meta.nChannels,1);
fields = {'nucLevel','nucLevel','nucLevel','nucLevel'};
for ci = 1:meta.nChannels
    im = imadjust(max(img{ci},[],3));
    cdata = newCellData.(fields{ci})(:,ci);
    ax(ci) = subplot_tight(2,2,ci);
    cellDataOverlay(im, XY, cdata, 0.02);
end
linkaxes(ax)
pause
close all

%% Process all positions and time points
nchannels = meta.nChannels;
nz = meta.nZslices;
opts.nuclearSegmentation = [];

tic
positions(meta.nPositions) = Colony();
for pidx = 1:meta.nPositions
    
    positions(pidx) = Colony(meta, pidx);
    positions(pidx).setRadius(350, meta.xres);
    cidx = ceil(pidx/meta.posPerCondition);
    positions(pidx).well = cidx;
    positions(pidx).dataChannels = 1:positions(pidx).nChannels;
    positions(pidx).filenameFormat = fnameformat;
    %iterate over all time points for each position
    for ti = 1:meta.nTime
        fprintf('t = %d\n',ti)
        %initialize cells for images + segmentations
        seg = cell(1,nchannels);
        img = cell(1,nchannels);
        %load images in each channel
        for cii = 1:nchannels
            ci = opts.dataChannels(cii);
            img{ci+1} = positions(pidx).loadImage(imdir,ci,ti);
        end
        %read nuclear segmentation
        name = fullfile(segdir,sprintf(bare,pidx-1,nucChannel,ti-1));
        segname = [name,'_FinalSegmentation.tif'];
        nucseg = false(size(img{1},1),size(img{1},2),nz);
        for zi = 1:nz
            nucseg(:,:,zi) = imread(segname,zi);
        end
        seg{nucChannel+1} = nucseg;
        
        %read cytoplasmic segmentation
        if ~isempty(opts.fgChannel)
            name = fullfile(segdir,sprintf(bare,pidx-1,opts.fgChannel,ti-1));
            segname = [name,'_Simple Segmentation.h5'];
            seg{opts.fgChannel+1} = ilastikRead(segname);
        end
        
        cellData = struct;
        for zi = 1:nz
            ims = cellfun(@(x) x(:,:,zi),img,'UniformOutput',false);
            segs = cell(1,nchannels);
            segs{nucChannel+1} = seg{nucChannel+1}(:,:,zi);
            if ~isempty(opts.fgChannel)
                segs{opts.fgChannel+1} = seg{opts.fgChannel+1}(:,:,zi);
            end
            opts.nuclearSegmentation = segs{nucChannel+1};
            [dbInfo, stats] = extractDataSingleFrame(ims,segs,opts);
            fields = fieldnames(stats);
            for fi = 1:length(fields)
                cellData(zi).(fields{fi}) = stats.(fields{fi});
            end
        end
        [newCellData, chain] = combineZsegmentation(cellData, meta, zopts);
        fields = fieldnames(newCellData);
        for fi = 1:length(fields)
            positions(pidx).cellData(ti).(fields{fi}) = newCellData.(fields{fi});
        end
    end
    
    positions(pidx).setCenter();
    positions(pidx).makeRadialAvgSeg();
    
    for ti = 1:meta.nTime
        positions(pidx).ncells(ti) = size(positions(pidx).cellData(ti).XY,1);
    end
    
    save(fullfile(dataDir,'positions.mat'),'positions')
end
toc

%% Or load positions if already processed
load(fullfile(dataDir,'positions.mat'))








