
clear all; close all;
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);

% raw data location, modify if not the same as the location of this script
dataDir = scriptPath; 
cd(dataDir);

%If there are multiple files for one well, the script picks the first file
%for each well for the overview, assuming that it is representative

% make an array of all the filenames
% remover FusionStitcher to get all images
listing = dir('*FusionStitcher.ims');

% do this if you haven't named your files like the well names
fileNames = {listing.name};
% to get the same structure as with regular makeOverview:
for i = 1:numel(fileNames)
    fileNames{i} = {fileNames{i}};
end

%Channel labels are given in the order that they were imaged
channelLabels = {{'DAPI','AP2y','SOX17','FOXA2'}};

%%
%Channel color picks which channels get which colors; if the first cell
%array of channelLabels is {'CDX2','DAPI','PDX1'} and the first array in
%channelColor is [3 1 2], then PDX1 is red, CDX2 is green, and DAPI is
%blue; i.e. channelOrders{i} = [indexForRed,indexForGreen,indexForBlue]
% 0 is gray, so [0 1 2] would make PDX1 gray
channelColor = {[0 0 0 0]};

%channelLabelIdx lets you pick which wells were imaged with which channels
%(or [] for no channels)
channelLabelIdx = num2cell(ones(size(fileNames)));

% channels starting from 1
channelsUsed = [1 2 3 4];

meta = struct('dataDir',dataDir,...
                        'channels', channelsUsed,...
                        'channelColor', {channelColor},...
                        'conditions',{fileNames},... % for this unstructured case
                        'fileNames',{fileNames},...
                        'channelLabels',{channelLabels},...
                        'channelLabelIdx',{channelLabelIdx});

%%
% read all the images
images = cell(fileNames);
for i = 1:size(fileNames,1)
    for j = 1:size(fileNames,2)
        for k = 1:numel(fileNames{i,j})
            images{i,j}{k} = readStack(fullfile(dataDir,fileNames{i,j}{k}));
        end
    end
end

%%
% determine intensity limits
tol = [0.03 0.03 0.03 0.03];
IlimAll = getIntensityLimits(images, meta, tol);

%%
% save previews for each channel separately

for i = 1:numel(channelsUsed)
    meta.channels = channelsUsed(i);
    postfix = ['_' channelLabels{1}{meta.channels}];
    savePreviews(images, meta, IlimAll, postfix)
end

%%
% make overlay of channels specified in channelsUsed

meta.channels = channelsUsed;
postfix = ['_' channelLabels{1}{meta.channels}];
savePreviews(images, meta, IlimAll, postfix)

    
