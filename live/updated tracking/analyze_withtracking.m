clear; close all;

% raw data location, modify if not the same as the location of this script
dataDir = 'Y:\Seth\200319_BMPlive_18h_63x_tracking';
mipDir = fullfile(dataDir,'MIP');
fnameformat = 'stitched_p%.4d_w%.4d_t%.4d.tif';

meta = Metadata(dataDir);

%%
% manual metadata
%-------------------

%since I am only looking at one stitched image right now
meta.nPositions = 2;
meta.nWells = 2;
meta.posPerCondition = 1;
meta.nTime = 50;

if isempty(meta.nPositions) || meta.nPositions ~= meta.nWells*meta.posPerCondition
    warning('total position doesnt match nWells x posPerCondition');
end

% in order of imaging (clockwise, which is not the layout of the dish)
meta.conditions = { 'Low Nuc Marker', 'High Nuc Marker' };

% set this to true if making an '8-well' loop through a 4-well
meta.loop4well = false;

% name your channels here if you want
meta.nucChannel = 0;
nucChannel = meta.nucChannel;
smadChannel = 1;
meta.channelLabel = {'Nuc','SMAD4'};

% specify if timelapse:
treatmentTime = 5; % frame of treatment
meta.timeInterval = '5 min';

%% set options to extract cell data
% notes:
% Andor channel counting starts at 0
% with cytoplasm: 'cytoplasmicLevels',true, 'cytoSize', 8 (example)

opts = struct(...
                    'dataChannels',     [nucChannel smadChannel],...
                    'fgChannel',        smadChannel,...
                    'cytoplasmicLevels',true,... 
                    'segmentationDir',  fullfile(dataDir,'MIP'),...
                    'nucShrinkage',     2,...
                    'bgMargin',         10,...
                    'segFG',            1,...
                    'saveCleanNucMask', true,...
                    'MIPidxDir',        fullfile(dataDir,'MIP'));

opts.cleanupOptions = struct('separateFused', true,...
    'clearBorder',true, 'minAreaStd', 1, 'minSolidity',0.95, 'minArea', 250);

%% check that the options are set right
close all
pidx = 2;
P = Position(meta);
P.filenameFormat = fnameformat;
P.setID(pidx);
segdir = fullfile(dataDir,'MIP');

ti = 1;
seg = P.loadSegmentation(segdir,nucChannel);
seg = seg(:,:,ti);
seg = seg == 1;
img = P.loadImage(dataDir,nucChannel,ti);
img = im2double(imadjust(max(img,[],3)));
newseg = nuclearCleanup(seg,opts.cleanupOptions);
seg = 0.5*double(seg); newseg = 0.5*double(newseg); img = 0.5*img;

imshow(cat(3,img+seg,img+newseg,img+newseg))


%% Do nuclear cleanup on each frame and save as tiffs

%still need to put this in

%% Load clean segmentations if previous step has already been completed

segDir = fullfile(dataDir, "Clean segmentations");
listing = dir(segDir);
names = {listing.name};
names = names(cellfun(@(x) contains(x,sprintf('_p%.4d',0)), names));
names = fullfile(segDir,names);
segs = cell(1,1,meta.nTime);
for ti = 1:meta.nTime
    segs{ti} = imread(names{ti});
end
segs = cell2mat(segs);

%% run the analysis on all positions and time points

tic
positions(meta.nPositions) = Position();
opts.tMax = meta.nTime;

segDir = fullfile(dataDir, "Clean segmentations");
listing = dir(segDir);

for pidx = 1:meta.nPositions
    disp(strcat("Position #",num2str(pidx)))
    names = {listing.name};
    names = names(cellfun(@(x) contains(x,sprintf('_p%.4d',pidx-1)),names));
    names = fullfile(segDir,names);
    segs = cell(1,1,meta.nTime);
    disp("Reading segmentations")
    for ti = 1:meta.nTime
        disp('.')
        segs{ti} = imread(names{ti}) == 1;
    end
    segs = cell2mat(segs);
    opts.nuclearSegmentation = segs;
    positions(pidx) = Position(meta);
    positions(pidx).filenameFormat = 'stitched_p%.4d_w%.4d_t%.4d.tif';
    positions(pidx).setID(pidx);
    positions(pidx).extractData(dataDir, opts);
%     positions(pi).makeAvgTimeTraces();
    save(fullfile(dataDir,'positions'), 'positions');
end
toc

%% load results if above block was run previously

load(fullfile(dataDir,'positions'));

%% load object classification and add to cellData
close all
listing = dir(mipDir);
names = {listing.name};
names = names(cellfun(@(x) contains(x,"Object Predictions"),names));
names = fullfile(mipDir,names);

for pidx = 1:meta.nPositions
    fnames = ...
        names(cellfun(@(x) contains(x,sprintf('_p%.4d',pidx-1)), names));
    with_time = ...
        fnames(cellfun(@(x) contains(x,sprintf('_t%.4d',0)), fnames));
    if ~isempty(with_time)
        seg = cell(1,1,meta.nTime);
        for ti = 1:meta.nTime
            fname =...
                fnames(cellfun(@(x) contains(x,sprintf('_t%.4d',ti-1)), fnames));
            fname = fname{1};
            seg{ti} = squeeze(h5read(fname,'/exported_data'))';
        end
        seg = cell2mat(seg);
    else
        fname = fnames{1};
        seg = h5read(fname, '/exported_data');
        seg = permute(squeeze(seg),[2,1,3]);
    end
%     for ti = 1:50
%         clf
%         imshow(double(seg(:,:,ti))/2);
%         pause 
%     end
    positions(pidx).addCellLabels(seg);
end

%% Check that labels were added correctly
close all
f = figure;
XY = {positions(2).cellData.XY};
labels = {positions(2).cellData.labels};
for ti = 1:meta.nTime
    clf
    img = positions(2).loadImage(dataDir,0,ti);
    img = imadjust(max(img,[],3));
    imshow(img)
    hold on
    scatter(XY{ti}(labels{ti} == 1,1),XY{ti}(labels{ti} == 1,2),36,'g','filled');
    scatter(XY{ti}(labels{ti} == 2,1),XY{ti}(labels{ti} == 2,2),36,'b','filled');
    scatter(XY{ti}(labels{ti} == 3,1),XY{ti}(labels{ti} == 3,2),36,'r','filled');
    title(num2str(ti))
    f.WindowState = 'maximized';
    pause
end

%% Set tracking parameters
%Parameters, constants
%max frame-frame cell movement distance in pixels; farther cells are ignored
max_linking_distance = 30; %~25-30 for 40x data; ~15 for 20x data
max_gap_distance = 70; %approximately double the max linking distance
max_gap_time = 5;

trackopts = struct(...
    'max_linking_distance', 50,...
    'max_gap_distance',     100,...
    'max_gap_time',         5,...
    'validate_links',       false,...
    'dataDir',              dataDir);

%% Do tracking for all time points
for pidx = 1:meta.nPositions
    try
        img = positions(pidx).loadImage(dataDir,nucChannel,1);
        trackopts.imsize = size(img);
        imsize = trackopts.imsize;
        disp(strcat("Image size = ",num2str(imsize(1)),"x",num2str(imsize(2))))
    catch
        try
            img = positions(pidx).loadImage(fullfile(dataDir,'MIP'),nucChannel,1);
            trackopts.imsize = size(img);
            imsize = trackopts.imsize;
            disp(strcat("Image size = ",num2str(imsize(1)),"x",num2str(imsize(2))))
        catch
            disp('Could not load image, using default size of 1024x1024')
            trackopts.imsize = [1024,1024];
        end
    end
    trees = create_tracks_v2(positions(pidx), trackopts);
    positions(pidx).tree = trees;
    save(fullfile(dataDir,'positions'), 'positions');
end

%% Visualize tracking results
close all

visopts = struct('tail_length', 20,...
    'write_to_tif', false,...
    'write_to_avi', true,...
    'default_tif', true,...
    'use_tif', true,...
    'use_seg', false,...
    'default_seg', true,...
    'default_name', true,...
    'visible', 'off',...
    'dataDir', dataDir);
%Making track visualization
tic
visualizeTracks5(positions(2), visopts);
toc






