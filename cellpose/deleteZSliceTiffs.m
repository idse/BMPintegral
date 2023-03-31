clear; close all; clc;

scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
dataDir = scriptPath;

%find subfolders of datadir with names ending with '_zslices'
folderlist = dir(fullfile(dataDir,'*_zslices'));

%iterate over folders and delete z slice and MIP tiff files and
%segmentation overlay images, but not segmentations (h5 and png files)
for ii = 1:length(folderlist)
    folder = fullfile(dataDir,folderlist(ii).name);
    %find tiff files corresponding to z slices
    imlist = dir(fullfile(folder,'*.tif'));
    %in case segmentation files are .tif instead of .png:
    imlist = imlist(cellfun(@(x) ~contains(x,'FinalSegmentation'),{imlist.name}));
    for fi = 1:length(imlist)
        name = fullfile(folder,imlist(fi).name);
        delete(name)
    end
    %find and delete segmentation overlay images
    imlist = dir(fullfile(folder,'*_SegOverlay*'));
    for fi = 1:length(imlist)
        name = fullfile(folder,imlist(fi).name);
        delete(name)
    end
end



