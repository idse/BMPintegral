clear; close all; clc;

scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
dataDir = scriptPath;

%find subfolders of datadir with names ending with '_zslices'
folderlist = dir(fullfile(dataDir,'*_zslices'));
%set nucChannel (in case it is ever not 0)
nucChannel = 0;
%identifier for nuclear channel images
nucid = sprintf('_w%.4d',nucChannel);

%iterate over folders and delete tiff files inside that include '_z' but
%not '_w0000' or 'MIP'
for ii = 1:length(folderlist)
    folder = fullfile(dataDir,folderlist(ii).name);
    %find tiff files corresponding to z slices
    imlist = dir(fullfile(folder,'*_z*.tif'));
    %exclude nuclear channel images
    imlist = imlist(cellfun(@(x) ~contains(x,nucid),{imlist.name}));
    %exclude MIP (I don't think this is strictly necessary since MIPs don't
    %include '_z', but just to be safe)
    imlist = imlist(cellfun(@(x) ~contains(x,'MIP'),{imlist.name}));
    for fi = 1:length(imlist)
        name = fullfile(folder,imlist(fi).name);
        delete(name)
    end
    %find an delete segmentation overlay images
    imlist = dir(fullfile(folder,'*_SegOverlay*'));
    for fi = 1:length(imlist)
        name = fullfile(folder,imlist(fi).name);
        delete(name)
    end
end



