clear; close all;

scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);

% raw data location, modify if not the same as the location of this script
dataDir = scriptPath; 

filelist = dir(fullfile(dataDir,'*FusionStitcher*.ims')); 

%convex decomposition options -> these are extremely permissive options
%tau1 is the threshold relative to cut length
tau1 = 0.7;
%tau2 is the absolute concavity threshold
tau2 = 1;
opts = struct('tau3',3);

%area threshold for getting rid of junk
minarea = 250;

nucChannel = 0;

for fi = 1:numel(filelist)
    
    fname = filelist(fi).name;
    s = strsplit(fname,{'_FusionStitcher','.ims'});
    prefix = [s{:}];
    
    tifname = strrep(strrep(fname,'FusionStitcher','Stitched'),'.ims','.tif');
    if exist(fullfile(dataDir,tifname),'file') == 2
        fname = tifname;
    end
    disp(fname)
    
    r = bfGetReader(fullfile(dataDir, fname));    
    nc = r.getSizeC();
    extension = fname(end-3:end);
    if strcmp(extension,'.tif')
        nz = r.getSizeZ();
    elseif strcmp(extension,'.ims')
        % -1 is only because fusion stitcher attaches a black z-slice at
        % the end; may have to be modified later
        nz = r.getSizeZ() - 1;
    end
    
    subDir = fullfile(dataDir,[prefix '_zslices']);
    pattern = [prefix '_w%.4d_z%.4d'];
    tic
    %figure
    %iterate over the z-slices
    for zi = 1:nz

        name = fullfile(subDir,sprintf(pattern,nucChannel,zi-1));
        fname = [name,'.tif'];
        segname = [name,'_Simple Segmentation.h5'];
        cpname = [name,'_cp_masks.png'];

        % read image and segmentations
        img = imread(fname);
        I = imadjust(img,stitchedlim(img,0.01));

        cpmask = imread(cpname);
        ilastikmask = ilastikRead(segname);
        mask = ilastikmask & ~(cpmask > 0);
        seg = imopen(mask,strel('disk',4));

        % do convex decomposition on remainder
        newseg = separate_fused(seg, tau1, tau2, opts);

        % combine the two masks
        L = uint16(labelmatrix(bwconncomp(newseg)));
        L(L>0) = L(L>0) + max(cpmask,[],'all');

        finalseg = bitor(cpmask,L);
        props = regionprops(finalseg,'PixelIdxList','Area');
        idxs = cell2mat({props([props.Area] < minarea).PixelIdxList}');
        finalseg(idxs) = 0;

        % save the resulting final mask
        savename = [name,'_FinalSegmentation.png'];
        imwrite(finalseg,savename);
        
        LRGB = label2rgb(finalseg,'hsv','k','shuffle');
        overlay = uint8(mat2gray(I)*255) + uint8(LRGB*0.5);
        savename = [name,'_SegOverlay.jpg'];
        imwrite(overlay,savename);%,'BitDepth',16)
% 
%         %visualization
%         im3 = visualize_nuclei_v2(finalseg,I);
%         clf
%         ax = gobjects(1,2);
%         ax(1) = subplot(1,2,1);
%         ax(1).Position([1,3]) = [0.01,0.485];
%         imshow(I)
%         ax(2) = subplot(1,2,2);
%         ax(2).Position([1,3]) = [0.505,0.48];
%         imshow(im3)
%         linkaxes(ax)
%         sgtitle(sprintf('%d',zi))
%         pause
    end
    toc
end










