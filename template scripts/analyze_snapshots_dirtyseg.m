clear all; close all;
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);

% raw data location, modify if not the same as the location of this script
%dataDir = scriptPath; 
dataDir = '/Users/idse/data_tmp/190513 pmax GFP lipo transfection test/Stitcher';

cd(dataDir);

% conditions table, follow layout of the plate
% this script assumes filenames start file well location, e.g. A1, B3,..
conditions = {'FW 1ul d1 75k','FW 1ul d1 125k','FW 1ul d1 175k';...
              'FW 1ul d2 50k','FW 1ul d2 83k','FW 1ul d2 100k';...
              'x','x','x';...
              'FW 2ul d2 50k','FW 2ul d2 83k','FW 2ul d2 100k';...
              'RV 1ul d1 75k','RV 1ul d1 125k','RV 1ul d1 175k'};

% make an array of all the filenames
listing = dir('*.ims');
[fileNames,wellNames] = getPlateFileNames(listing, conditions);

%%

efficiency = [];
confluence = [];

for i = 1:size(fileNames,1)
    for j = 1:size(fileNames,2)
        for k = 1:numel(fileNames{i,j})
        
            fileName = fileNames{i,j}{k};

            if exist(fileName)

                nuclearSegFname = [fileName(1:end-4) '_MIP_w0000_Simple Segmentation.h5'];
                GFPSegFname = [fileName(1:end-4) '_MIP_w0001_Simple Segmentation.h5'];

                image = readStack(fullfile(dataDir,fileName));
                adjustedNuclei = imadjust(mat2gray(image(:,:,1)));
                adjustedGFP = imadjust(mat2gray(image(:,:,2)));

                nuclearSeg = h5read(fullfile(dataDir,'MIP',nuclearSegFname),'/exported_data');
                nuclearSeg = squeeze(nuclearSeg)'==2;

                GFPSeg = h5read(fullfile(dataDir,'MIP',GFPSegFname),'/exported_data');
                GFPSeg = squeeze(GFPSeg)'==2;

                GFP = mat2gray(image(:,:,2));
                adjustedGFP = imadjust(GFP,stretchlim(GFP,[0.1 0.95]));
                GFPSegOutline = GFPSeg - imerode(GFPSeg, strel('disk',2));
                overlay = cat(3,adjustedNuclei + GFPSegOutline,adjustedGFP + GFPSegOutline,nuclearSeg);
                imshow(overlay)

                imwrite(overlay(1:1024,1:1024,:),fullfile(dataDir,[fileName(1:end-4) '_seg.tif']));

                efficiency(i,j,k) = sum(GFPSeg & nuclearSeg)/sum(nuclearSeg);
                confluence(i,j,k) = sum(sum(imclose(nuclearSeg,strel('disk',8))))/numel(nuclearSeg);
            end
        end
    end
end

%%

for i = 1:size(efficiency,1)
    for j = 1:size(efficiency,2)
        disp([wellNames{i,j} ' ' conditions{i,j} '  :  '  num2str(mean(efficiency(i,j,:))) '  :  ' num2str(mean(confluence(i,j,:)))]);
    end
end

            