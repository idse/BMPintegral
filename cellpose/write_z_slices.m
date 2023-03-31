clear; close all; clc;

scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);

% raw data location, modify if not the same as the location of this script
dataDir = scriptPath; 
% set nucChannel; it is almost always 0
nucChannel = 0;
%find FusionStitcher files
filelist = dir(fullfile(dataDir,'**/*FusionStitcher*.ims'));
%don't process stuff in 'sorted' folder
filelist = filelist(...
    cellfun(@(x) ~contains(x,[filesep,'sorted']),{filelist.folder}));
%don't process files in any subfolder named 'temp'
filelist = filelist(...
    cellfun(@(x) ~contains(x,[filesep,'temp']),{filelist.folder}));
%create a binary vector to keep track of which files have stitching errors
badColonyMask = false(numel(filelist),1);
for fi = 1:numel(filelist)
    
    fname = filelist(fi).name;
    s = strsplit(fname,{'_FusionStitcher','.ims'});
    prefix = [s{:}];
    
    subDir = filelist(fi).folder;
    writeDir = fullfile(subDir,[prefix '_zslices']);
    mkdir(writeDir)

    outfnamePattern = [prefix '_w%.4d_z%.4d.tif'];
    MIPpattern = [prefix,'_MIP_w%.4d.tif'];

    disp(['writing z slices for ' fname]);
    tic
    r = bfGetReader(fullfile(subDir, fname));
    r.setSeries(0);
    
    nZslices = r.getSizeZ();
    nChannels = r.getSizeC();
    MIPs = cell(nChannels,1);
%     ci = nucChannel;
    for ci = 0:nChannels-1
        for zi = 0:nZslices-1
            im = bfGetPlane(r, r.getIndex(zi,ci,0)+1);
            if isempty(MIPs{ci+1})
                MIPs{ci+1} = im;
            else
                MIPs{ci+1} = max(MIPs{ci+1},im);
            end
            if ci==nucChannel
                writename = fullfile(writeDir,sprintf(outfnamePattern,ci,zi));
                percentBlank = 100*sum(im==0,'all')/numel(im);
                if percentBlank > 90 && zi == nZslices-1
                    %this is likely because the last z slice is completely
                    %black, just don't use that one
                    break
                elseif percentBlank > 10
                    %indicates an issue with the Fusion stitcher
                    warnMessage =...
                        strcat("Fusion stitching error detected,",...
                        " using lab stitching for ", prefix);
                    warning(warnMessage)
                    badColonyMask(fi) = true;
                    break
                end
                imwrite(im,writename)
            end
        end
        if badColonyMask(fi)
            break
        end
        MIPname = fullfile(writeDir,sprintf(MIPpattern,ci));
        imwrite(MIPs{ci+1},MIPname)
    end
    toc
end

badidxs = find(badColonyMask);

disp('Colonies with stitching errors:')
disp({filelist(badidxs).name}')

%% Do our stitching for colonies with detected errors
montageGridSize = [3 3];
montageOverlap = 20;
%assumes 1024x1024 images
pixelOverlap = round(1024*montageOverlap/100);
conditionPositionFusion = [7 2 9; 4 5 6; 1 8 3];

for ind = 1:numel(badidxs)
    fi = badidxs(ind);
    fname = filelist(fi).name;
    s = strsplit(fname,{'_FusionStitcher','.ims'});
    prefix = [s{:}];
    subDir = filelist(fi).folder;
    writeDir = fullfile(subDir,[prefix '_zslices']);
    if ~isfolder(writeDir)
        mkdir(writeDir)
    else
        %delete any existing contents first
        rmdir(writeDir,'s')
        mkdir(writeDir)
    end
    outfnamePattern = [prefix '_w%.4d_z%.4d.tif'];
    
    %find images
    if isempty(s{2})
        %if fusion stitcher has no _F%d identifier, find filenames in the
        %main dataDir
        imdir = subDir;
        list = dir(fullfile(imdir,[prefix,'*.ims']));
        %exclude the FusionStitcher
        list = list(cellfun(@(x) ~contains(x,'FusionStitcher'),{list.name}));
    else
        %if fname is of the form prefix_FusionStitcher_F%d.ims, look for
        %individual images in a subfolder
        idxstr = strsplit(s{2},'_F');
        %should there be an error if there is not _F in the string?
        idx = str2double(idxstr{end});
        imdir = fullfile(subDir,[s{1},sprintf('_Field_%d',idx+1)]);
        list = dir(fullfile(imdir,'*.ims'));
    end
    %read all images in each grid position
    nucims = cell(montageGridSize);
    imgs = cell(montageGridSize);
    if numel(list) ~= numel(conditionPositionFusion)
        error('Expected to find %d individual position IMS files, but found %d instead\n',...
            numel(conditionPositionFusion),numel(list))
    end
    for pidx = 1:numel(conditionPositionFusion)
        readname = fullfile(imdir,list(pidx).name);
        [jj,ii] = ind2sub(montageGridSize, conditionPositionFusion(pidx));
        imgs{ii,jj} = readStack(readname);
        nucims{ii,jj} = max(squeeze(imgs{ii,jj}(:,:,nucChannel+1,:)),[],3);
    end
    %determine # of channels, # of z slices from stacks of images
    nZslices = size(imgs{1,1},4);
    nChannels = size(imgs{1,1},3);
    %determine stitching offsets based on nuclear channel MIPs
    maxPixelOverlap = round(size(nucims{1,1}, 2)/4);
    offset = cross_corr(nucims{1,1},nucims{1,2}, 2, maxPixelOverlap);
    if abs(offset/pixelOverlap - 1) > 0.25
        warning('detected pixel overlap very different from nominal value, sticking with nominal value');
    else
        pixelOverlap = offset;
    end
    upperleft = registerImageGrid_v3(nucims, pixelOverlap);
    %stitch the nuclear images to determine the image size
    [stitched, ~] =...
        stitchImageGridIntensityWeighted(upperleft, nucims);
    xSize = size(stitched,1); ySize = size(stitched,2);
    clear nucims %mostly to save ram
    
    %iterate over channels and z stacks to do stitching and write stitched
    %stack as a 4D tiff
    stitched = zeros(xSize,ySize,nChannels,nZslices,class(imgs{1}));
    for ci = 1:nChannels
        for zi = 1:nZslices
            ims = cellfun(@(x) x(:,:,ci,zi),imgs,'UniformOutput',false);
            [stitch, ~] =...
                stitchImageGridIntensityWeighted(upperleft, ims);
%             [stitch, ~] =...
%                 stitchImageGridWeightedAvg(upperleft, ims);
            stitched(:,:,ci,zi) = stitch;
        end
    end
    %write as a tif with 'FusionStitcher' replaced by 'Stitched'
    tifname = strrep(strrep(fname,'FusionStitcher','Stitched'),'.ims','.tif');
    tifname = fullfile(subDir,tifname);
    writeTiffStack(stitched,tifname)
    %write nuclear z slices
    ci = nucChannel;
    for zi = 1:nZslices
        im = stitched(:,:,ci+1,zi);
        writename = fullfile(writeDir,sprintf(outfnamePattern,ci,zi-1));
        imwrite(im,writename)
    end
    %write MIPs
    for ci = 1:nChannels
        MIPname = fullfile(writeDir,sprintf(MIPpattern,ci-1));
        MIP = max(squeeze(stitched(:,:,ci,:)),[],3);
        imwrite(MIP,MIPname)
    end
    
end
