clear all; close all;

scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);

% raw data location, modify if not the same as the location of this script
dataDir = scriptPath; 

% load metadata that was specified in makeMIP
load(fullfile(dataDir,'meta.mat'));

% manual metadata
%-------------------

meta.nPositions = meta.nWells*meta.posPerCondition;
meta.nucChannel = 0;
meta.channelLabel = {'H2B','pERK','BRA','ERK-KTR'};
meta.filenameFormat = 'stitched_p%.4d_w%.4d_t%.4d.tif';

save(fullfile(dataDir,'meta.mat'),'meta');

nucChannel = 0;
pERKChannel = 1;
BRAChannel = 2;
ERKKTRChannel = 3;

%% load an image and segmentation, display as sanity check

ci = 0;
ti = 1;
pi = 1;

P = Colony(meta);
P.setID(pi);
im = P.loadImage(dataDir, ci, ti);
seg = P.loadSegmentation(fullfile(dataDir,'MIP'), ci);

MIP = max(im,[],3);
imshow(imadjust(MIP),[]);

figure, 
imshow(seg(:,:,ti),[])

%% set options for extracting nuclear and cytoplasmic levels

% externally I will have all indices starting at 1
% the Andor offset to start at 0 will be internal

opts = struct(...
                    'dataChannels',     0:3,... %'fgChannel',        S4Channel,...
                    'cytoplasmicLevels',true,... 
                    'segmentationDir',  fullfile(dataDir,'MIP'),...
                    'tMax',             1,...
                    'nucShrinkage',     2,...
                    'cytoSize',         5,...
                    'cytoMargin',       0,...
                    'bgMargin',         10,...
                    'segFG',            1,...
                    'MIPidxDir',        fullfile(dataDir,'MIP'),...
                    'NCRcutoff',        []);

% opts.cleanupOptions = struct('separateFused', true,...
%     'clearBorder',true, 'minAreaStd', 0, 'minSolidity',0.95, 'minArea',100);

% TODO:
% add explanation for what dirtyOptions does
opts.dirtyOptions = struct('s', 5, 'areacutoff', 100, 'toobigNstd', 3); % try toobigNstd = 3

%% check that the standard options are set right

pi = 1;
P = Colony(meta, pi);
ti = 1;
opts.tMax = ti;
dbInfo = P.extractData(dataDir, opts);

P.cellDataOverview(ti);

im = P.loadImage(dataDir, nucChannel, ti);
MIP = max(im,[],3);

figure,
A = 0.5*imadjust(mat2gray(MIP));
s = 0.3;
imshow(cat(3, A + dbInfo.nucmaskraw*s, A + s*dbInfo.nucmask, A + s*dbInfo.cytmask));

opts.tMax = meta.nTime;

%% show segmentation overlay another way

figure,
segoverlay = nuclearSegmentationOverlay(MIP, dbInfo.nucmask);
imshow(segoverlay);

%% set colony parameters

radiusMicron = 350;
P.setRadius(radiusMicron, meta.xres);

%% make radial average and plot

options = struct(   'nucChannels', [pERKChannel, BRAChannel, nucChannel],...
                    'cytChannels', pERKChannel);
plotRadialProfiles(P, meta, options);

%% make a normalized radial plot

options = struct(   'nucChannels', [pERKChannel, BRAChannel, nucChannel],...
                    'cytChannels', pERKChannel,...
                    'normalize', true);
plotRadialProfiles(P, meta, options);

%% for instructional purposes without the function

P.makeRadialAvgSeg();

nuc_profile = P.radialProfile.NucAvgSeg;
cyt_profile = P.radialProfile.CytAvgSeg;
r = P.radialProfile.BinEdges(1:end-1)*meta.xres;
%r = (r(1:end-1)+r(2:end))/2;

%figure,
clf
hold on
plot(r, nuc_profile(:,pERKChannel+1),'LineWidth',2, 'Color', 'g')
plot(r, cyt_profile(:,pERKChannel+1),'LineWidth',2, 'Color', 'b')
plot(r, nuc_profile(:,BRAChannel+1),'LineWidth',2, 'Color', 'r')
plot(r, nuc_profile(:,nucChannel+1),'LineWidth',2, 'Color', 'k')
hold off
legend({'nuc ppERK','cyt ppERK', 'BRA', 'H2B'});

% make it pretty
fgc = 'k';
bgc = 'w';
graphbgc = 1*[1 1 1]; 
graphfgc = 'r';
fs = 24;

xlabel(['radius ( ... )'], 'FontSize',fs,'FontWeight','Bold','Color',fgc)
ylabel('intensity', 'FontSize',fs,'FontWeight','Bold','Color',fgc);
    
set(gcf,'color',bgc);
set(gca, 'LineWidth', 2);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')
set(gca,'XColor',fgc);
set(gca,'YColor',fgc);
set(gca,'Color',graphbgc);

%% normalize radial plots

P.makeRadialAvgSeg();

nuc_profile = P.radialProfile.NucAvgSeg;
cyt_profile = P.radialProfile.CytAvgSeg;

for ci = 1:size(nuc_profile,2)
    nuc_profile(:,ci) = (nuc_profile(:,ci) - min(nuc_profile(:,ci)))./(max(nuc_profile(:,ci)) - min(nuc_profile(:,ci)));
    cyt_profile(:,ci) = (cyt_profile(:,ci) - min(cyt_profile(:,ci)))./(max(cyt_profile(:,ci)) - min(cyt_profile(:,ci)));
end

r = P.radialProfile.BinEdges(1:end-1)*meta.xres;
%r = (r(1:end-1)+r(2:end))/2;

%figure,
plot(r, nuc_profile(:,pERKChannel+1),'LineWidth',2, 'Color', 'g')
hold on
plot(r, cyt_profile(:,pERKChannel+1),'LineWidth',2, 'Color', 'b')
plot(r, nuc_profile(:,BRAChannel+1),'LineWidth',2, 'Color', 'r')
plot(r, nuc_profile(:,nucChannel+1),'LineWidth',2, 'Color', 'k')
hold off
legend({'nuc ppERK','cyt ppERK', 'BRA', 'H2B'}, 'Location','NorthWest');

% make it pretty
fgc = 'k';
bgc = 'w';
graphbgc = 1*[1 1 1]; 
graphfgc = 'r';
fs = 24;

xlabel(['radius ( ... )'], 'FontSize',fs,'FontWeight','Bold','Color',fgc)
ylabel('intensity', 'FontSize',fs,'FontWeight','Bold','Color',fgc);
    
set(gcf,'color',bgc);
set(gca, 'LineWidth', 2);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')
set(gca,'XColor',fgc);
set(gca,'YColor',fgc);
set(gca,'Color',graphbgc);


%% run the analysis on all colonies

opts.tMax = meta.nTime;

tic
positions(meta.nPositions) = Colony();

for pi = 1:meta.nPositions

    positions(pi) = Colony(meta);
    positions(pi).setID(pi);
    
    positions(pi).setRadius(350, meta.xres);
    
    dbInfo = positions(pi).extractData(dataDir, opts);
    
    positions(pi).makeRadialAvgSeg();
    save(fullfile(dataDir,'positions'), 'positions');
    
    % save segmentation overlay
    time = 1;
    im = positions(pi).loadImage(dataDir, meta.nucChannel, time);
	MIP = max(im,[],3);
    segoverlay = nuclearSegmentationOverlay(MIP, dbInfo.nucmask);
    imwrite(segoverlay, fullfile(dataDir, sprintf('segOverlay_p%.4d.png', pi-1)));
    
    % save radial profile
    options = struct(   'nucChannels', [pERKChannel, BRAChannel, nucChannel],...
                        'cytChannels', pERKChannel,...
                        'normalize', false);
    plotRadialProfiles(positions(pi), meta, options);
    saveas(gcf, fullfile(dataDir, ['plot_' num2str(pi) '.png']));
    
    % save radial profile normalized
    options.normalize = true;
    plotRadialProfiles(positions(pi), meta, options);
    saveas(gcf, fullfile(dataDir, ['plotNorm_' num2str(pi) '.png']));
end
toc

%% make average radial plot for all colonies

options = struct(   'nucChannels', [pERKChannel, BRAChannel, nucChannel],...
                    'cytChannels', pERKChannel,...
                    'normalize', false,...
                    'std', true);
plotRadialProfiles(positions, meta, options);

%%
options.normalize = true;
plotRadialProfiles(positions, meta, options);

%% compare variation in colonies 


%% scatter plots colored for radius


%% radial profiles without segmenting individual nuclei for comparison


%% colonies over time

