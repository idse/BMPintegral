clear; close all; clc

%% consolidate cp masks
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
dataDir = scriptPath;
% bare = 'stitched_p%.4d_w0000_t0000_z%.4d_cp_masks.png';
bare = 'stitched_p%.4d_w%.4d_t%.4d';
savebare = 'stitched_p%.4d_w%.4d_t%.4d_cp_masks.tif';
nucChannel = 0;

listing = dir(fullfile(dataDir,['stitched_p*',sprintf('_w%.4d_t%.4d.tif',nucChannel,0)]));
npos = length(listing);
listing = dir(fullfile(dataDir,[sprintf('stitched_p%.4d_w%.4d',0,nucChannel),'_t*.tif']));
listing = listing(~cellfun(@(x) contains(x,'FinalSegmentation') | contains(x,'cp_masks'),{listing.name}));
ntime = length(listing);

%% reformat cellpose masks
outputdir = fullfile(dataDir, 'cp_masks');
if ~exist(outputdir,'dir'), mkdir(outputdir); end

for ii = 1:npos
    prefix = sprintf(bare,ii-1,nucChannel,0);
    listing = dir(fullfile(dataDir,[prefix,'*cp_masks.png']));
    nz = length(listing);
    disp(prefix)
    fprintf('nz = %d\n',nz)
    for ti = 1:ntime
        prefix = sprintf(bare,ii-1,nucChannel,ti-1);
        savename = fullfile(dataDir,sprintf(savebare,ii-1,nucChannel,ti-1));
        for zi = 1:nz
            fprintf('.')
            if zi == 1
                mode = 'overwrite';
            else
                mode = 'append';
            end
            name = strcat(prefix,'_z000',num2str(zi-1),'_cp_masks.png');%[prefix,sprintf('_z%.4d_cp_masks.png',zi-1)];
            readname = fullfile(dataDir,name);
            img = imread(readname);
            imwrite(img,savename,'WriteMode',mode)

            writename = fullfile(outputdir,[prefix,sprintf('_z%.4d_cp_masks.png',zi-1)]);
            movefile(readname,writename)
        end
        fprintf('\n')
    end
end

%% combine ilastik, cellpose
bare = 'stitched_p%.4d_w%.4d_t%.4d';
%convex decomposition parameters
tau1 = 0.7;
tau2 = 1;
opts = struct('tau3',3);
%area threshold for getting rid of junk
minarea = 50; maxarea = 800;

overlaydir = fullfile(dataDir, 'SegOverlays');
if ~exist(overlaydir,'dir'), mkdir(overlaydir); end
tic
for ii = 1:npos
    for ti = 1:ntime
        fname = sprintf(bare,ii-1,nucChannel,ti-1);
        fname = fullfile(dataDir,fname);
        cpname = [fname,'_cp_masks.tif'];
        ilastikname = [fname '_Simple Segmentation.h5'];
        imname = [fname,'.tif'];
        ilastikseg = ilastikRead(ilastikname);
        nz = size(ilastikseg,3);
        sname = [fname,'_FinalSegmentation.tif'];
        for zi = 1:nz
            if zi == 1
                mode = 'overwrite';
            else
                mode = 'append';
            end
            img = imread(imname,zi);
            I = imadjust(img,stitchedlim(img,0.01));
            cpmask = imread(cpname,zi);
            ilastikmask = ilastikseg(:,:,zi);
            %difference of the masks
            mask = ilastikmask & ~(cpmask > 0);
            seg = imopen(mask,strel('disk',4));
            % do convex decomposition on remainder
            newseg = separate_fused(seg, tau1, tau2, opts);
            % combine the two masks
            L = uint16(labelmatrix(bwconncomp(newseg)));
            L(L>0) = L(L>0) + max(cpmask,[],'all');

            finalseg = bitor(cpmask,L);
            finalseg = finalseg > 0 & ~boundarymask(finalseg);
            props = regionprops(finalseg,'PixelIdxList','Area');
            idxs = cell2mat({props([props.Area] < minarea | [props.Area] > maxarea).PixelIdxList}');
            finalseg(idxs) = 0;

            finalseg = uint8(finalseg);
            
            imwrite(finalseg,sname,'WriteMode',mode)
            
            [~,name,~] = fileparts(fname);
            overlay = visualize_nuclei_v2(finalseg>0,I);
            savename = [name, sprintf('_z%.4d',zi-1), '_SegOverlay.jpg'];
            savename = fullfile(overlaydir,savename);
            imwrite(overlay,savename);
        end
    end
end
toc



