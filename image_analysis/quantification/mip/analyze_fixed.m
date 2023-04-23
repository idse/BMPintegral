clear; close all; clc

%% setup
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
dataDir = scriptPath;
mipDir = fullfile(dataDir,'MIP');
segdir = mipDir;

listing = dir(fullfile(mipDir,'*_MIP_p*_w0000.tif'));
names = {listing.name};
prefixes = cell(size(names));
for ii = 1:length(prefixes)
    I = strfind(names{ii},'_MIP_');
    prefixes{ii} = names{ii}(1:I-1);
end
prefixes = unique(prefixes);
np = NaN(size(prefixes));
for ii = 1:length(prefixes)
    np(ii) = sum(contains(names,prefixes{ii}));
end
ids = [0,cumsum(np)];
didxs = NaN(length(names),1); dpis = NaN(length(names),1);
for ii = 1:length(np)
    didxs(ids(ii)+1:ids(ii+1)) = ii;
    dpis(ids(ii)+1:ids(ii+1)) = 1:np(ii);
end

conditions = {'mTeSR,0hr','mTeSR,0hr',...
    'LDN0,6hr','LDN10,6hr','LDN30,6hr',...
    'LDN0,12hr','LDN10,12hr','LDN30,12hr','LDN100,12hr','LDN300,12hr','mTeSR,12hr',...
    'LDN0,24hr','LDN10,24hr','LDN30,24hr','LDN100,24hr','LDN300,24hr','mTeSR,24hr',...
    'LDN0,36hr','LDN10,36hr','LDN30,36hr','LDN100,36hr','LDN300,36hr','mTeSR,36hr',...
    'LDN0,48hr','LDN10,48hr','LDN30,48hr','LDN100,48hr','LDN300,48hr','mTeSR,48hr',...
    };

npos = length(names);
ncond = length(conditions);
ppc = npos/ncond;
channelLabel = {'DAPI','ISL1','SOX2','NANOG'};
nc = length(channelLabel);
nucChannel = 0;

manualMeta = struct();
manualMeta.plateLayoutNames = false;
manualMeta.loop4well = false;
%conditions
manualMeta.conditions = conditions;
% number of experimental conditions
manualMeta.nWells = ncond;
manualMeta.posPerCondition = ppc;
manualMeta.nPositions = npos;

manualMeta.nucChannel = nucChannel;
manualMeta.channelLabel = channelLabel;
% complete the rest of the metadata automatically and make metadata object
meta = Metadata(dataDir, manualMeta);

save(fullfile(dataDir,'meta'), 'meta');


opts = struct(...
                    'cytoplasmicLevels',true,...
                    'cytoSize',         4,...
                    'cytoMargin',       4,...
                    'nucShrinkage',     2,...
                    'bgMargin',         0,...
                    'nucChannel',       nucChannel,...
                    'junkmask',         []);

%% run quantification on all positions

opts.suppressOutput = false;
positions(npos) = Position();
tic
for pidx = 1:npos
    %which dataset/directory is this from, and what is its index among
    %images from that subset
    didx = didxs(pidx); dpi = dpis(pidx);
    %what condition, and what index within that condition?
    cidx = ceil(pidx/ppc); cpi = pidx - (cidx - 1)*ppc;
    prefix = prefixes{didx};
    
    segname = sprintf([prefix,'_MIP_p%.4d_w%.4d_FinalSegmentation.tif'],dpi-1,nucChannel);
    objectname = sprintf([prefix,'_MIP_p%.4d_w%.4d_Object Predictions.h5'],dpi-1,nucChannel);
    fname = sprintf([prefix,'_F%.2d.ims'],dpi-1);
    
    seg = imread(fullfile(segdir,segname)) == 1;
    objectseg = ilastikRead(fullfile(segdir,objectname));
    masks = makeNucCytMasks(seg,[],opts);
    
    fprintf('Position %d of %d\n',pidx-1,npos-1)
    positions(pidx) = Position(meta, dpi);
    positions(pidx).filenameFormat = [prefix,'_F%.2d.ims'];
    
    imgs = cell(1,nc);
    for ci = 1:nc
        imgs{ci} = positions(pidx).loadImage(dataDir,ci-1,1);
    end
    %make masks, read out nuclear and cytoplasmic levels
    celldata = readoutFromMasks(imgs,masks,opts);
    positions(pidx).cellData = celldata;
    positions(pidx).ncells = size(celldata.XY,1);
    positions(pidx).addCellLabels(objectseg);
    
    save(fullfile(dataDir,'positions.mat'),'positions')
end
toc


%% load data if already processed
load(fullfile(dataDir,'positions.mat'),'positions')



