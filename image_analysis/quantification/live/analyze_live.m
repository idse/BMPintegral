clear; close all;
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);

% raw data location, modify if not the same as the location of this script
dataDir = scriptPath;

% metadata
%-------------------

manualMeta = struct();

% experimental conditions
manualMeta.conditions = {'30k,50ng/mL'};
manualMeta.nWells = length(manualMeta.conditions);
manualMeta.posPerCondition = 6;
manualMeta.nPositions = manualMeta.nWells*manualMeta.posPerCondition;

% set this to true for scenario 2: if names have A1, A2 etc to encode plate layout:
manualMeta.plateLayoutNames = false;

% set this to true if making an '8-well' loop through a 4-well
manualMeta.loop4well = false;

manualMeta.nucChannel = 0;
manualMeta.channelLabel = {'H2B','SMAD4'};

% this should be automatic but is not in metadata due to fusion bug
manualMeta.timeInterval = '10 min';

% complete the rest of the metadata automatically and make metadata object
meta = Metadata(dataDir, manualMeta);

save(fullfile(dataDir,'meta'), 'meta');

% additional parameters
%-------------------

% name your channels here if you want
nucChannel = meta.nucChannel;
smadChannel = 1;

fnameformat = 'stitched_p%.4d_w%.4d_t%.4d.tif';


%% extract nuclear and cytoplasmic levels
opts = struct(...
                    'dataChannels',     0:1,...
                    'cytoplasmicLevels',true,...
                    'fgChannel',        smadChannel,...
                    'cytoSize',         5,...
                    'segmentationDir',  fullfile(dataDir,'MIP'),...
                    'MIPidxDir',        fullfile(dataDir,'MIP'),...
                    'nucShrinkage',     2,...
                    'bgMargin',         6,...
                    'segFG',            1); % Ilastik class defining the foreground

%nuclear segmentation cleanup options
opts.cleanupOptions = struct('separateFused', false,...
    'clearBorder',true,...
    'minAreaStd', 1,...
    'minSolidity',0.95,...
    'minArea', 400,...
    'openSize', 1,...
    'fillholes', true);

%nuclear splitting via convex decomposition
tau1 = 0.7; %relative concavity threshold
tau2 = 1.5; %absolute concavity threshold
minArea = 600; %don't cut objects smaller than this
maxArea = 3000; %don't consider objects larger than this to be dividing nuclei

opts.decompopts = struct(...
    'flag',         false,...
    'tau3',         6,...
    'useMinArea',   true,...
    'minArea',      minArea,...
    'tau1',         tau1,...
    'tau2',         tau2,...
    'maxArea',      maxArea,...
    'ignoreholes',  false);

fgChannels = [1, 2, 3];
divChannels = [2, 3];


%% check the nuclear cleanup options
close all
pidx = 1;
times = 1;
segdir = opts.segmentationDir;
list = dir([segdir,'/*Object Predictions.h5']);
names = fullfile(segdir,{list.name});

segname = names(cellfun(@(x) contains(x,sprintf('p%.4d',pidx-1)), names));
allseg = cell(1,1,length(times));
for ti = 1:length(times)
    allseg{ti} = ilastikRead(segname{times(ti)});
end
allseg = cell2mat(allseg);

f = figure('WindowState','maximized');
for ti = 1:length(times)
    clf
    seg = allseg(:,:,ti);
    tic
    newseg = newNuclearCleanup(seg, fgChannels, divChannels, opts);
    toc
    
    P = Position(meta);
    P.filenameFormat = fnameformat;
    P.setID(pidx);
    img = max(P.loadImage(dataDir,nucChannel,times(ti)),[],3);
    img = im2double(imadjust(img,stitchedlim(img)));
    nimg = visualize_nuclei(newseg,img);
    
    ax2 = subplot(1,2,2);
    imshow(nimg)
    title('final segmentation with image')
    
    ax2.Position([1,3]) = [0.51, 0.48];
    
    ax1 = subplot(1,2,1);
    imshowpair(ismember(seg,fgChannels),newseg)
    title('changes from cleanup')
    ax1.Position([1,3]) = [0.01, 0.48];
    linkaxes([ax1,ax2])
    sgtitle(sprintf('Time = %d', times(ti)))
    pause
end
close(f)

%% Test options for extractData (run above block first)
close all
tempopts = opts;
tempopts.nuclearSegmentation = newseg;
tempopts.tMax = ti;
P.nTime = ti;
dbInfo = extractData(P, dataDir, tempopts);
P.cellDataOverview(ti);

im = mat2gray(max(P.loadImage(dataDir, opts.fgChannel, ti),[],3));
im = 0.7*imadjust(im,stitchedlim(im));

figure
% A = 0.5*imadjust(mat2gray(MIP));
s = 0.3;
A = img;
ax1 = subplot_tight(1,3,1);
imshow(cat(3, A + dbInfo.nucmaskraw*s, A + s*dbInfo.nucmask, A + s*dbInfo.cytmask));
ax2 = subplot_tight(1,3,2);
imshow(cat(3, im + dbInfo.nucmaskraw*s, im + s*dbInfo.nucmask, im + s*dbInfo.cytmask));
ax3 = subplot_tight(1,3,3);
imshow(im/0.7);

linkaxes([ax1,ax2,ax3])


opts.tMax = meta.nTime;

imsize = size(A);

%% Do nuclear cleanup on each frame and save as tiffs

writeDir = fullfile(dataDir,'Clean segmentations');
if ~exist(writeDir,'dir'), mkdir(writeDir); end
pattern = 'clean_seg_p%.4d_w%.4d.tif';

segdir = opts.segmentationDir;
list = dir([segdir,'/*Object Predictions.h5']);
names = fullfile(segdir,{list.name});

tic
for pidx = 1:meta.nPositions
    newsegs = zeros([imsize(1),imsize(2),meta.nTime],'uint8');
    disp(strcat("Processing position ", num2str(pidx)))
    %process, clean up segmentation at each frame, and convert to uint8
    parfor ti = 1:meta.nTime
        segname = names{cellfun(@(x) contains(x,sprintf('p%.4d_w%.4d_t%.4d',pidx-1,nucChannel,ti-1)),names)}; %#ok<PFBNS>
        seg = ilastikRead(segname);    
        newseg = newNuclearCleanup(seg, fgChannels, divChannels, opts);
        newseg = uint8(newseg);
        newsegs(:,:,ti) = newseg;
    end
    %write files
    for ti = 1:meta.nTime
        if mod(ti,round(meta.nTime/20)) == 0
            fprintf('.')
        end
        writename = fullfile(writeDir,sprintf(pattern, pidx-1, nucChannel));
        if ti == 1
            imwrite(newsegs(:,:,ti), writename);
        else
            imwrite(newsegs(:,:,ti), writename, 'WriteMode', 'append')
        end
    end
    fprintf('\n')
end
toc

%% run the analysis on all positions and time points
close all
positions(meta.nPositions) = Position();
segdir = opts.segmentationDir;

tic
segDir = fullfile(dataDir, "Clean segmentations");
listing = dir(fullfile(segDir,'/*.tif'));
names = {listing.name};
opts.tMax = meta.nTime;
% for all conditions (aka wells)
for ci = 1:meta.nWells
    
    % for all positions belonging to that condition
    for cpi = 1:meta.posPerCondition
        
        pi = meta.posPerCondition*(ci-1) + cpi;
        disp(strcat("Position #", num2str(pi)))
        segnames = names(cellfun(@(x) contains(x,sprintf('_p%.4d',pi-1)),names));
        segnames = fullfile(segDir,segnames);
        segs = cell(1,1,meta.nTime);
        disp("Reading segmentations")
        for ti = 1:opts.tMax
            if mod(ti,round(opts.tMax/20)) == 0
                fprintf('.')
            end
            segs{ti} = imread(segnames{1},ti) == 1;
        end
        fprintf('\n')
        segs = cell2mat(segs);
        opts.nuclearSegmentation = segs;
        
        positions(pi) = Position(meta, pi);
        positions(pi).filenameFormat = fnameformat;
        
        debugInfo = positions(pi).extractData(dataDir, opts);
        positions(pi).makeAvgTimeTraces();
        save(fullfile(dataDir,'positions'), 'positions');
    end
end
toc

%% load results if above block was run previously

load(fullfile(dataDir,'positions'));

%% load object classification and add to cellData
close all
%set segDir to the full path of the folder containing object classification
%results
segDir = fullfile(dataDir, "MIP");
listing = dir(segDir);
names = {listing.name};
names = names(cellfun(@(x) contains(x,"Object Predictions"),names));
names = fullfile(segDir,names);

tic
for pidx = 1:meta.nPositions
    fprintf('Position #%d\n',pidx)
    fnames = ...
    names(cellfun(@(x) contains(x,sprintf('_p%.4d',pidx-1)), names));
    seg = cell(1,1,meta.nTime);
    for ti = 1:meta.nTime
        seg{ti} = ilastikRead(fnames{ti});
    end
    seg = cell2mat(seg);
    positions(pidx).addCellLabels(seg);
end
fprintf('\n')
toc

save(fullfile(dataDir,'positions'), 'positions');

%% Check that labels were added correctly
close all

pidx = 1;
XY = {positions(pidx).cellData.XY};
labels = {positions(pidx).cellData.labels};
for ti = [1 25 75 150 250]
    img = max(positions(pidx).loadImage(dataDir,nucChannel,ti),[],3);
    img = imadjust(img,stitchedlim(img));
    cla
    imshow(img)%,'Parent',ha);
    hold on
    scatter(XY{ti}(labels{ti} == 1,1),XY{ti}(labels{ti} == 1,2),36,'g','filled');
    scatter(XY{ti}(labels{ti} == 2,1),XY{ti}(labels{ti} == 2,2),36,'b','filled');
    scatter(XY{ti}(labels{ti} == 3,1),XY{ti}(labels{ti} == 3,2),36,'r','filled');
    title(num2str(ti))
%     f.WindowState = 'maximized';
    pause
end


