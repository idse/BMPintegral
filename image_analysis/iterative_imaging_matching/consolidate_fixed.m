%% Load data, make LineageTrace object
clear; close all;
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);

% raw data location, modify if not the same as the location of this script
dataDir = scriptPath;

baseDir = 'Y:\Seth\220712_Smad4GFP_BMP_IWP2_MP_live';
liveDir = fullfile(baseDir,'fixed_round1');
fixedDir = fullfile(baseDir,'fixed_round2');

liveMeta = load(fullfile(liveDir,'meta.mat'));
liveMeta = liveMeta.meta;

fixedMeta = load(fullfile(fixedDir,'meta.mat'));
fixedMeta = fixedMeta.meta;

%initialize lineage trace
lt = LineageTrace(liveDir,fixedDir);
npos = length(lt.live_position);

%get image size
A = lt.fixed_position(1).loadImage(fixedDir,0,1);
imsize = size(A,[1,2]);
clear A
%number of time points in live imaging
ntime = lt.live_position(1).nTime;
nucChannel = 0;

%% Map fixed round 1 to fixed round 1
close all
positionIdx = 1:npos;
maxDist = 6; %um
mapPoints(lt, maxDist, struct('alignmentMethod','automatic',...
    'positionIdx',positionIdx,'useZ',true));

for ii = positionIdx
    lt.checkAlignment(ii);
    cleanSubplot(18)
    title(sprintf('Position #%d',ii))
    pause
    clf
end
close all

% save(fullfile(dataDir,'lt.mat'),'lt')

%% Or load existing LineageTrace object
load(fullfile(dataDir,'lt.mat'))

%% add cell data from second round to first round
load(fullfile(liveDir,'positions.mat'),'positions')
load(fullfile(liveDir,'meta.mat'),'meta')

channelLabels = [meta.channelLabel, fixedMeta.channelLabel];
c1 = meta.nChannels; c2 = fixedMeta.nChannels;
nchannels = c1 + c2;

meta.channelLabel = channelLabels;
meta.nChannels = nchannels;

fields = {'cytLevel','nucLevel','NCratio'};

for ii = 1:npos
    ncells = positions(ii).ncells;
    mapped_idxs = lt.mapped_idxs{ii};
    fixedIdxs = find(mapped_idxs > 0);
    liveIdxs = mapped_idxs(fixedIdxs);
    for fi = 1:length(fields)
        field = fields{fi};
        A = NaN(ncells,nchannels);
        A1 = lt.live_position(ii).cellData.(field);
        if ~isempty(A1)
            A(:,1:c1) = A1;
        end
        A2 = lt.fixed_position(ii).cellData.(field);
        if ~isempty(A2)
            A(liveIdxs,c1+1:nchannels) = A2(fixedIdxs,:);
        end
        positions(ii).cellData.(field) = A;
    end
    bg = [lt.live_position(ii).cellData.background,...
        lt.fixed_position(ii).cellData.background];
    positions(ii).cellData.background = bg;
end

%%
save(fullfile(liveDir,'combined_positions.mat'),'positions')
save(fullfile(liveDir,'combined_meta.mat'),'meta')

