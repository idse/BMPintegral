clear; close all; clc

%%
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
dataDir = scriptPath;

listing = dir(fullfile(dataDir,'*_MIP_*_w0000.tif'));

npos = length(listing);
ntime = 1;

%convex decomposition parameters
tau1 = 0.8;
tau2 = 1.2;
opts = struct('tau3',3,'minArea',1500,'useMinArea',false);
%area threshold for getting rid of junk
minarea = 70;
maxarea = 2500;

%% combine ilastik, cellpose

overlaydir = fullfile(dataDir, 'SegOverlays');
if ~exist(overlaydir,'dir'), mkdir(overlaydir); end

tic
for ii = 1:npos
    for ti = 1:ntime
%         fname = fullfile(dataDir,sprintf(bare,ii-1,0));
        fname = fullfile(dataDir,listing(ii).name(1:end-4));
        cpname = [fname,'_cp_masks.png'];
        ilastikname = [fname '_Simple Segmentation.h5'];
        imname = [fname,'.tif'];
        ilastikmask = ilastikRead(ilastikname);
        sname = [fname,'_FinalSegmentation.tif'];
        
        [filepath,name,ext] = fileparts(imname);
        disp(name)
        
        img = imread(imname);
        cpmask = imread(cpname);
        mask = ilastikmask & ~(cpmask > 0);
        
        seg = imopen(mask,strel('disk',4));
        % do convex decomposition on remainder
        newseg = separate_fused(seg, tau1, tau2, opts);
        % combine the two masks
        L = uint16(labelmatrix(bwconncomp(newseg)));
        L(L>0) = L(L>0) + max(cpmask,[],'all');

        finalseg = bitor(cpmask,L);
        finalseg = bwmorph(finalseg > 0 & ~boundarymask(finalseg),'thicken',1);
        props = regionprops(finalseg,'PixelIdxList','Area');
        idxs = cell2mat({props([props.Area] < minarea | [props.Area] > maxarea).PixelIdxList}');
        finalseg(idxs) = 0;
        
        finalseg = uint8(finalseg);
        
        % save the resulting final mask
        imwrite(finalseg,sname)
        if ti == 1 || ti == ntime
            I = imadjust(img,stitchedlim(img,0.01));
            overlay = visualize_nuclei_v2(finalseg>0,I);
            savename = [name, '_SegOverlay.jpg'];
            savename = fullfile(overlaydir,savename);
            imwrite(overlay,savename);
        end
        
    end
end

toc