clear; close all; 

scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);

% raw data location, modify if not the same as the location of this script
dataDir = scriptPath; 
cd(dataDir);

manMeta = struct();
manMeta.nucChannel = 0;
manMeta.channelLabel = {'DAPI','ISL1','EOMES','SOX17'};

manMeta.conditions = {'ctrl'};
manMeta.posPerCondition = 5;
manMeta.nWells = numel(manMeta.conditions);
manMeta.nPositions = manMeta.posPerCondition*manMeta.nWells;

filelist = {};
filelist{1} = dir(fullfile(dataDir,'**/210330_42_HR_ISL1-EOMES-SOX17*FusionStitcher.ims'));

posfile = fullfile(dataDir,['positions' [manMeta.channelLabel{:} '.mat']]);
if exist(posfile)
	load(posfile, 'positions');
    disp(['loading ' posfile]);
end

statsfile = fullfile(dataDir,['stats' [manMeta.channelLabel{:} '.mat']]);
if exist(statsfile)
	load(statsfile, 'stats');
    disp(['loading ' statsfile]);
end

countfile = fullfile(dataDir,['counts' [manMeta.channelLabel{:} '.mat']]);
if exist(countfile)
	load(countfile, 'counts');
    disp(['loading ' countfile]);
end

metafile = fullfile(dataDir,['meta'  [manMeta.channelLabel{:} '.mat']]);
if exist(metafile)
    load(metafile,'meta');
    disp(['loading ' metafile]);
else
    fname = filelist{1}(1).name;
    subdir = filelist{1}(1).folder;
    meta = Metadata(fullfile(subdir,fname), manMeta);
    save(metafile, 'meta');
end

%% extractData
opts = struct('cytoSize',4,'cytoMargin',2,'IoU',0.5);

positions(manMeta.nPositions) = Colony();

% for all conditions (aka wells)
tic
for ci = 1:manMeta.nWells
    
    % for all positions belonging to that condition
    for cpi = 1:manMeta.posPerCondition
    
        pi = manMeta.posPerCondition*(ci-1) + cpi;
        
        fname = filelist{ci}(cpi).name;
        subdir = filelist{ci}(cpi).folder;
        
        meta = Metadata(fullfile(subdir,fname), manMeta);
        % -1 is only because fusion stitcher attaches a black z-slice at the end
        % may have to be modified later
        meta.nZslices = meta.nZslices - 1;
        %if there is a stitched 4D tiff with a name matching this one, use
        %that instead of the FusionStitcher file
        tifname = strrep(strrep(fname,'FusionStitcher','Stitched'),'.ims','.tif');
        if exist(fullfile(subdir,tifname),'file') == 2
            fname = tifname;
            r = bfGetReader(fullfile(subdir,tifname));
            %update metadata
            meta.nZslices = r.getSizeZ(); %don't subtract 1 for stitched tiff stacks
            meta.xSize = r.getSizeX();
            meta.ySize = r.getSizeY();
            fnameformat = meta.filenameFormat;
            meta.filenameFormat =...
                strrep(strrep(fnameformat,'FusionStitcher','Stitched'),'.ims','.tif');
        end
        
        % id ~= pi in this case
        if ~isempty(regexp(fname, '_F[0-9]','Once'))
            extension = fname(end-3:end);
            s = strsplit(fname,{'_F',extension});
            id = uint16(str2double(s{end-1})) + 1;
        else
            id = pi;
        end
        fprintf('Position %d\n',pi)
        disp(fname)
        positions(pi) = Colony(meta, id);
        positions(pi).setRadius(350, meta.xres);
        positions(pi).well = ci;
        positions(pi).dataChannels = 1:positions(pi).nChannels;
        
        positions(pi).cellData = extractDataMultiZ(subdir, fname, meta, opts);
        positions(pi).setCenter();
        positions(pi).makeRadialAvgSeg();
        positions(pi).ncells = size(positions(pi).cellData.XY,1);
    end
end
toc
save(posfile, 'positions');

%% make xsections of segmentation
%save 1 of these per condition
tic
for ci = 1:manMeta.nWells
    cpi = 1;
    pi = manMeta.posPerCondition*(ci-1) + cpi;
    fname = positions(pi).filename;
    prefix = fname(1:end-4);
    extension = fname(end-3:end);
    img = readStack(fullfile(dataDir,fname));
    
    LMfname = fullfile(dataDir,[prefix '_LabelMatrix3D.mat']);
    LM = load(LMfname);
    LM = LM.LM;
    
    s = strsplit(fname, {'_Stitched','_FusionStitcher',extension});
    meta = load(fullfile(dataDir,[[s{:}] '_zslices'],'meta.mat'));
    meta = meta.meta;
    nuclearim = squeeze(img(:,:,nucChannel+1,1:meta.nZslices));
    %change yinds and xinds if preferred
    yinds = 1400;
    xinds = 1:size(img,2);
    makeCrossSection(nuclearim,LM,meta,yinds,xinds)
    saveas(gcf, fullfile(dataDir, [prefix  'segmentationXsection_y' num2str(yinds) '.png'])); 
    close all
end
toc

%% make intensity radial profile

% find shared limits for all positions
for pi = 1:numel(positions)
    if pi == 1
        maxvals = max(positions(pi).radialProfile.NucAvgSeg(1:end-10,:));
        minvals = min(positions(pi).radialProfile.NucAvgSeg(1:end-10,:));
    elseif pi < 5
        maxvals = max(max(positions(pi).radialProfile.NucAvgSeg(1:end-10,:)), maxvals);
        minvals = min(min(positions(pi).radialProfile.NucAvgSeg(1:end-10,:)), minvals);
    end
end
nuclimits = [minvals' maxvals'];
        
figure('Position',[0 0 500 500])
options = struct(   'nucChannels', 1:3,...
                    'normalize', true, 'std', true, 'FontSize', 30, 'legend',false);%, ...
                    %'nuclimits', nuclimits);

for ci = 1:manMeta.nWells
    
    pi = manMeta.posPerCondition*(ci-1)+1:manMeta.posPerCondition*ci;
    clf
    plotRadialProfiles(positions(pi), meta, options);
    %ylim([0 4000])
    ylim([0 1.1]);
    axis square

    saveas(gcf, fullfile(dataDir, [meta.conditions{ci} '_' meta.channelLabel{2:4} '_radialProfiles_v2.png'])); 
end

%% make distributions out of processed data

stats = cellStats(positions, meta, positions(1).dataChannels(2:end));
tolerance = 0.01;
nbins = 50;
stats.makeHistograms(nbins, tolerance);
confidence = 0.95;
whichthreshold = [1 1 2];
conditions = 1:numel(meta.conditions);
stats.getThresholds(confidence, conditions);%, whichthreshold);
%stats.exportCSV(meta); export to CSV
save(statsfile, 'stats');

clf
conditionIdx = 1:numel(meta.conditions);
for channelIdx = 1:4
    options = struct('channelIndex',channelIdx, 'cytoplasmic', false,...
                                'cumulative',false, 'time', 1,...
                            'conditionIdx',conditionIdx,...
                            'titlestr',meta.channelLabel{channelIdx}...
                            );
    clf                
    stats.plotDistributionComparison(options)
    saveas(gcf, fullfile(dataDir, [meta.channelLabel{channelIdx(1)} '_dist.png'])); 
end
%ylim([0 0.02])
%xlim([0 1000])

%%
condi = 1;
i = 2;
j = 3;

% SAMPLE QC PLOT
X = stats.nucLevel{condi}(:,i);
Y = stats.nucLevel{condi}(:,j);
X = log(1+X/mean(X));
Y = log(1+Y/mean(Y));
scatter(X, Y, 2, stats.sample{1})


%% scatter plot subpopulation on image
figure, 

normIdx = 1;
pi = 1;
ci = 3;

condi = ceil(pi / meta.posPerCondition);
cpi = pi - (condi-1)*meta.posPerCondition;
subDir = filelist{condi}(cpi).folder;

channelThresholds = stats.thresholds;

s = strsplit(positions(pi).filename,{'_FusionStitcher','.ims'});
prefix = [s{:}];

zdir = [prefix '_zslices'];
img = imread(fullfile(subDir, zdir, sprintf([prefix '_MIP_w%.4d.tif'], ci-1)));
%img = max(positions(pi).loadImage(dataDir,ci-1,1),[],3);
img = imadjust(img,stitchedlim(img));
%img = mat2gray(img, [150 4000]);

[X,Y] = meshgrid(1:size(img,2),1:size(img,1));
R = sqrt((X - positions(pi).center(1)).^2 + (Y - positions(pi).center(2)).^2);
mask = R > positions(pi).radiusPixel*1.05;
img(mask) = max(img(:));

nucLevel = positions(pi).cellData.nucLevel;
background = positions(pi).cellData.background;
Ncells = positions(pi).ncells;

positive = {};
for i = 1:meta.nChannels
    positive{i} = nucLevel(:,i) - background(i) > channelThresholds(i);
end

%figure,
imshow(img)
hold on
%scatter(positions(pi).cellData.XY(:,1),positions(pi).cellData.XY(:,2),'filled')
XY = positions(pi).cellData.XY(positive{ci},:);
scatter(XY(:,1),XY(:,2),20,'filled','r')
% XY = positions(pi).cellData.XY(AP2Cp,:);
% scatter(XY(:,1),XY(:,2),'x','g')
XY = positions(pi).cellData.XY(positive{2},:);
%scatter(XY(:,1),XY(:,2),'x','b')
scatter(positions(pi).center(1),positions(pi).center(2),1000,'.','g')
hold off
title(meta.channelLabel{ci})

for pi = 1:numel(positions)
    positions(pi).cellData.channelThresholds = channelThresholds;
end
save(fullfile(dataDir,'positions'), 'positions');


%% radial probability of being positive

Pall = {};
xall = {};
channelThresholds = [0 500 500 800];
stats.thresholds =  channelThresholds;%[0 500 400 700];
for condi = 1:numel(meta.conditions)
    combos = {};%[2 4],[2 3],[3 4]};
    [P,x,Pstd] = radialPositive(stats, positions, condi, meta,combos);
    Pall{condi} = P;
    Pstdall{condi} = Pstd;
    xall{condi} = x;
    ylim([0 1]);
    saveas(gcf, fullfile(dataDir, ['radialpositive_' meta.conditions{condi} '_', meta.channelLabel{2:end} '.png']));
    close;
end

%% count populations
            
combo = [3 2 4];
%combo = [4 2];
conditionsidx = 1:5;%6:9;
combosubset = [1 3];
counts = countPopulations(positions, meta, stats, dataDir, combo, conditionsidx);
save(countfile, 'counts');

%% make pretty scatter plot

close all;
options = struct();
options.conditionIdx = 1%:numel(meta.conditions);
options.channelCombos = {[3 2]};%, [4 3], [4 2]};
options.channelThresholds = stats.thresholds;
options.channelMax = exp([3 3 3.5 3])-1;
options.log1p = true;
options.radiusMicron = positions(1).radiusMicron;
options.conditionsCombined = false;

scatterMicropattern(stats, meta, dataDir, options)

%% pretty channel combo pie
   
options = struct();
options.channels = [2 3 4];
options.pieOrder = [3 1 2];
% tolerances in original order of channels
% manMeta.channelLabel = {'DAPI','AP2C','NANOG','PRDM1'};
%manMeta.channelLabel = {'DAPI','AP2C','EOMES','SOX17'};
options.tol = [0.01 0.99; 0.1 0.99; 0.8 0.99; 0.7 0.995];

posidx = [6 5 7 8]; % the first position sets the lookup table

for pi = posidx
    
    ci = ceil(pi/meta.posPerCondition);
    cpi = mod(pi-1, meta.posPerCondition)+1;
    fname = filelist{ci}(cpi).name;
    subDir = filelist{ci}(cpi).folder;

    if pi == posidx(1)
        Ilim = micropatternPieVis(subDir, positions(pi), options);
    else
        options.Ilim = Ilim;
        micropatternPieVis(subDir, positions(pi), options);
    end
end

%% multi condition pie
  
options = struct();
options.channels = [2 3 4];
% tolerances in original order of channels
% manMeta.channelLabel = {'DAPI','AP2C','NANOG','PRDM1'};
options.tol = [0.01 0.99; 0.1 0.99; 0.01 0.99; 0.7 0.995];
 
conditionIdx = [1 2 4 5 6];

options.positionIdx = (conditionIdx-1)*meta.posPerCondition + 1;

micropatternPieVisConditions(dataDir, positions, options, meta)






