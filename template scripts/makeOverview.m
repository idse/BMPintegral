
clear all; close all;
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);

% raw data location, modify if not the same as the location of this script
dataDir = scriptPath; 
cd(dataDir);

% conditions table, follow layout of the plate
% this script assumes filenames start file well location, e.g. A1, B3,..
% i.e. conditions{1,1} is the condition for filenames starting with A1,
% conditions{1,2} is for files with A2, conditions{2,1} is for files with
% B1, etc.
conditions = {'IGF 17nM', 'LY 10uM','RTKi 10uM', 'Pi103 50nM' 'SM';...
}; 

%If there are multiple files for one well, the script picks the first file
%for each well for the overview, assuming that it is representative

% make an array of all the filenames
listing = dir('*.ims');
[fileNames,wellNames] = getPlateFileNames(listing, conditions);

%Channel labels are given in the order that they were imaged
channelLabels = {{'DAPI','pAKT'}};

%%
%Channel color picks which channels get which colors; if the first cell
%array of channelLabels is {'CDX2','DAPI','PDX1'} and the first array in
%channelColor is [3 1 2], then PDX1 is red, CDX2 is green, and DAPI is
%blue; i.e. channelOrders{i} = [indexForRed,indexForGreen,indexForBlue]
% 0 is gray, so [0 1 2] would make PDX1 gray
channelColor = {[0 2]};

%channelLabelIdx lets you pick which wells were imaged with which channels
%(or [] for no channels)
channelLabelIdx = { 1, 1, 1, 1, 1;...
};

% channels starting from 1
channelsUsed = [1 2];

meta = struct('dataDir',dataDir,...
                        'channels', channelsUsed,...
                        'channelColor', {channelColor},...
                        'fileNames',{fileNames},...
                        'conditions',{conditions},...
                        'channelLabels',{channelLabels},...
                        'channelLabelIdx',{channelLabelIdx});

%%
% read all the images
images = cell(size(conditions));
for i = 1:size(conditions,1)
    for j = 1:size(conditions,2)
        for k = 1:numel(fileNames{i,j})
            images{i,j}{k} = readStack(fullfile(dataDir,fileNames{i,j}{k}));
            DAPIgray = true;
        end
    end
end

%%
% determine intensity limits
tol = 0.01;
IlimAll = getIntensityLimits(images, meta, tol);

%%
% save previews
savePreviews(images, meta, IlimAll)


