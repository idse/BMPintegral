clear all; close all;
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);

% raw data location, modify if not the same as the location of this script
dataDir = scriptPath;
cd(dataDir);

% metadata
%-------------------

manualMeta = struct();

% number of experimental conditions
manualMeta.nWells = 18;

% naming conventions:
% 
% the script needs to know how to organize files into conditions, there are
% two ways:
% 1) multiposition acquisition, ending in filenames X_F00, X_F01, X_FN
% in this case, way specify positionsPerCondition, and number of conditions
% and they should multiply to give N.
% 2) snapshots: name can be whatever, as long as you add a coordinate to
% indicate the condition, e.g. A1_X, A1_Y, A2_Z, or X_A1_Y, Z_A1_W. 
% These are labeled as well coordinates on a plate but don't
% have to correspond to an actual well coordinate. What is important is to
% start with A1 and have a contiguous block of wells. 
%
% set this to true for scenario 2: if names have A1, A2 etc to encode plate layout:
manualMeta.plateLayoutNames = true;

% set this to true if making an '8-well' loop through a 4-well
manualMeta.loop4well = false;

% as plate layout
manualMeta.conditions = {'0 min','10 min','30 min','60 min','120 min','240 min';...
                        '+RI 0 min','+RI 10 min','+RI 30 min','+RI 60 min','+RI 120 min','+RI 240 min';...
                        '+MEKi 30min','+MEKi 240min','C10A100RI (bra ctrl)','E6T 0min','E6T 10min','E6T 30min'};

manualMeta.nucChannel = 0;
manualMeta.channelLabel = {'DAPI','BRA','pERK'};

% this should be automatic but is not in metadata due to fusion bug
manualMeta.timeInterval = '10 min';

% complete the rest of the metadata automatically and make metadata object
meta = Metadata(dataDir, manualMeta);

save(fullfile(dataDir,'meta'), 'meta');

% additional parameters
%-------------------

% name your channels here if you want
nucChannel = meta.nucChannel;
pERKChannel = 2;

% specify if timelapse:
treatmentTime = 30; % frame of treatment

%% extract nuclear and cytoplasmic levels

% notes:
% Andor channel counting starts at 0
% with cytoplasm: 'cytoplasmicLevels',true, 'cytoSize', 8 (example)

% add this when using analysis with multiple z:
% 'MIPidxDir',        segDir,...

% add something like this if segmenting the cells
%'fgChannel',        pERKChannel,...

opts = struct(...
                    'dataChannels',     0:2,...
                    'cytoplasmicLevels',true,... 
                    'cytoSize',         4,...
                    'segmentationDir',  fullfile(dataDir,'MIP'),...
                    'nucShrinkage',     2,...
                    'bgMargin',         10,...
                    'segFG',            1); % Ilastik class defining the foreground

opts.cleanupOptions = struct('separateFused', true,...
    'clearBorder',true, 'minFusedArea', 1500, 'minSolidity',0.95,...
    'erodeSize',10,'minArea', 150);

%% check that the options are set right

pi = 1;
P = Position(meta, pi);
ti = 1;
opts.tMax = ti;
dbInfo = P.extractData(dataDir, opts);

P.cellDataOverview(ti);

im = P.loadImage(dataDir, pERKChannel, ti);
MIP = max(im,[],3);

figure,
A = 0.5*imadjust(mat2gray(MIP));
s = 0.3;
imshow(cat(3, A + dbInfo.nucmaskraw*s, A + s*dbInfo.nucmask, A + s*dbInfo.cytmask));

opts.tMax = meta.nTime;

%%
% (optional) try out the nuclear cleanup settings on some frame:
nucseg = dbInfo.segall{nucChannel+1};
cleannucseg = nuclearCleanup(nucseg(:,:,ti), opts.cleanupOptions);
segoverlay = nuclearSegmentationOverlay(MIP, cleannucseg);
imshow(segoverlay);

%% run the analysis on all positions and time points

tic
positions(meta.nPositions) = Position();

% for all conditions (aka wells)
for ci = 1:meta.nWells
    
    % for all positions belonging to that condition
    for cpi = 1:meta.posPerCondition

        pi = meta.posPerCondition*(ci-1) + cpi;
        positions(pi) = Position(meta, pi);

        dbInfo = positions(pi).extractData(dataDir, opts);
        positions(pi).makeAvgTimeTraces();
        save(fullfile(dataDir,'positions'), 'positions');

        % save segmentation overlay, this is just QC, can be commented out
        time = opts.tMax;
        im = positions(pi).loadImage(dataDir, meta.nucChannel, time);
        MIP = max(im,[],3);
        segoverlay = nuclearSegmentationOverlay(MIP, dbInfo.nucmask);
        imwrite(segoverlay, fullfile(dataDir, sprintf('segOverlay_p%.4d.png', pi-1)));

         % overlay some quantified property on an image, optional, can be
         % commented out
        pERKnuc = positions(pi).cellData(end).nucLevel(:,pERKChannel + 1); %index starts from 1, channel from 0
        pERKcyt = positions(pi).cellData(end).cytLevel(:,pERKChannel + 1); %index starts from 1, channel from 0
        im = max(positions(pi).loadImage(dataDir, pERKChannel, opts.tMax), [], 3);
        tol = 0.02;
%         cellDataOverlay(im, positions(pi).cellData(end).XY, pERKnuc, tol);
%         saveas(gcf, fullfile(dataDir, sprintf('perkNucValuesOverlay_p%.4d.png',pi)));
        cellDataOverlay(im, positions(pi).cellData(end).XY, pERKcyt, tol);
        saveas(gcf, fullfile(dataDir, sprintf('perkCytValuesOverlay_p%.4d.png',pi-1)));
    end
end
toc

%% load results if above block was run previously

load(fullfile(dataDir,'positions'));

% % to concatenate two files (if microscope had to be stopped)
% positions_before = load(fullfile(dataDir,'positions_before'));
% positions_before = positions_before.positions;
% positions = concatenatePositions(positions_before, positions);
% save(fullfile(dataDir,'positions'), 'positions');

%% make figure of the time traces (if timelapse data)

options = struct('wellnrs', 1:4, 'channel', pERKChannel,...
                'mode','N',...
                'yLabel', 'Bra:H2B::Cer',...
                'minNCells', 10,...
                'dataChannels', opts.dataChannels,...
                'treatmentTime', treatmentTime,... 
                'signalingRange', []);

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
legend(meta.conditions, 'Location','northwest');
export_fig(fullfile(dataDir,'timeTrace_multipleConditions.png'),'-native -m2');


%% show distributions in different conditions 

% make distributions out of processed data
stats = cellStats(positions, meta, opts.dataChannels);
tolerance = 0.03;
nbins = 50;
stats.makeHistograms(nbins, tolerance);

%%
% plot distribution and save figure
channelIndex = pERKChannel+1;

FGFtimeseriesIdx = [3 1 4 7 10 13 16 12 15];
FGFRItimeseriesIdx = [3 2 5 8 11 14 17 12 15];

for ti = 1

    for cum = [true false]
        clf
        % can leave time option out if static images
        options = struct('channelIndex',channelIndex, 'cytoplasmic', false,...
                        'cumulative',cum, 'time', ti,...
                        'conditionIdx',FGFtimeseriesIdx,...
                        'titlestr','24h C3+A100, then FGF2 -> pERK'...
                        );
        stats.plotDistributionComparison(options)
        ylim([0 0.3]);
        
        filename = ['distOverlay_' meta.channelLabel{channelIndex}];
        if options.cumulative
            filename = [filename 'cum'];
            ylim([0 1]);
        end
        %saveas(gcf, fullfile(dataDir, sprintf([filename '_t%.4d.png'],ti)));    
    end
end


%% make graphs for snapshot time series

% size(cytLevelMean) = size(conditions) x nChannels 

cytlevelMean = zeros([size(meta.conditions) numel(opts.dataChannels)]);
cytlevelStd = zeros([size(meta.conditions) numel(opts.dataChannels)]);

for ci = 1:meta.nWells

    [i,j] = ind2sub(size(meta.conditions),ci);

    cytlevelTot = [];

    for cpi = 1:meta.posPerCondition

        pi = meta.posPerCondition*(ci-1) + cpi;
        cytlevelTot = cat(1, cytlevelTot, positions(pi).cellData.cytLevel);
    end

    cytlevelMean(i,j,:) = nanmean(cytlevelTot);
    cytlevelStd(i,j,:) = nanstd(cytlevelTot);
end


%%

FGFtimeseriesIdx = [1 4 7 10 13 16];
FGFRItimeseriesIdx = [2 5 8 11 14 17];
MEKiIdx = [3 6];
E6TIdx = [12 15 18];

% time points
t = [0 10 30 60 120 240];

tlist = {t, t, [30 240], [0 10 30]}; 
idxlist = {FGFtimeseriesIdx, FGFRItimeseriesIdx, MEKiIdx, E6TIdx};



clf 
hold on
for ci = pERKChannel+1
    for i = 1:numel(idxlist)
        cavg = cytlevelMean(:,:,ci);
        cstd = cytlevelStd(:,:,ci);

        cavg = cavg(idxlist{i});
        cstd = cstd(idxlist{i});

        errorbar(tlist{i}, cavg, cstd, '-x', 'LineWidth',2); 
    end
end
hold off

fs = 16;
set(gca, 'LineWidth', 2);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')
xlabel('time (min)');
ylabel('ppERK (a.u.)');
xlim([min(t)-10, max(t) + 10]);
legend('C+A 24h; FGF2', 'C+A 24h; FGF2 RI', 'MEKi', 'E6T; FGF2');

saveas(gcf, 'timeTrace_snapshots.png');


%% variability between  images

% indices belonging to different conditions
FGFtimeseriesIdx = [1 4 7 10 13 16];
FGFRItimeseriesIdx = [2 5 8 11 14 17];
MEKiIdx = [3 6];
E6TIdx = [12 15 18];

% time points
t = [0 10 30 60 120 240];
tlist = {t, t, [30 240], [0 10 30]}; 
idxlist = {FGFtimeseriesIdx, FGFRItimeseriesIdx, MEKiIdx, E6TIdx};

% collect means
% size(cytLevelMean) = size(conditions) x posPerCondition x nChannels 
cytlevelMean = zeros([size(meta.conditions) meta.posPerCondition numel(opts.dataChannels)]);
cytlevelStd = zeros([size(meta.conditions) meta.posPerCondition numel(opts.dataChannels)]);

for ci = 1:meta.nWells

    [i,j] = ind2sub(size(meta.conditions),ci);

    cytlevelTot = [];

    for cpi = 1:meta.posPerCondition

        pi = meta.posPerCondition*(ci-1) + cpi;
        
        cytlevelMean(i,j,cpi,:) = nanmean(positions(pi).cellData.cytLevel);
        cytlevelStd(i,j,cpi, :) = nanstd(positions(pi).cellData.cytLevel);
    end

    
end

% make figure combined
clf 
colors = lines(4);
hold on
for ci = pERKChannel+1
for k = 1:meta.posPerCondition
    for i = 1:numel(idxlist)
        
        cavg = cytlevelMean(:,:,k,ci);
        cstd = cytlevelStd(:,:,k,ci);

        cavg = cavg(idxlist{i});
        cstd = cstd(idxlist{i});

        errorbar(tlist{i}, cavg, cstd, '-x', 'LineWidth',1,'Color', colors(i,:)); 
    end
end
end
hold off

fs = 16;
set(gca, 'LineWidth', 2);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')
xlabel('time (min)');
ylabel('ppERK (a.u.)');
xlim([min(t)-10, max(t) + 10]);
legendstr = {'C+A 30h; FGF2', 'C+A 30h; FGF2 RI', 'MEKi', 'E6T; FGF2'};
legend(legendstr);

saveas(gcf, 'timeTrace_snapshots_posvar.png');

%%
% make figure

for i = 1:numel(idxlist)
    
    clf 
    colors = lines(4);
    hold on
    
    for ci = pERKChannel+1
        for k = 1:meta.posPerCondition
        
            cavg = cytlevelMean(:,:,k,ci);
            cstd = cytlevelStd(:,:,k,ci);

            cavg = cavg(idxlist{i});
            cstd = cstd(idxlist{i});

            errorbar(tlist{i} + k, cavg, cstd, '-x', 'LineWidth',1,'Color', colors(i,:)); 
        end
    end
    fs = 16;
    set(gca, 'LineWidth', 2);
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
    xlabel('time (min)');
    ylabel('ppERK (a.u.)');
    xlim([min(t)-10, max(t) + 10]);
    ylim([500 4000]);
    legend('C+A 30h; FGF2', 'C+A 30h; FGF2 RI', 'MEKi', 'E6T; FGF2');

    saveas(gcf, ['timeTrace_' legendstr{i} '.png']);
    
    hold off
end






