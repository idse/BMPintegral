clear all; close all;
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);

% raw data location, modify if not the same as the location of this script
dataDir = scriptPath; 
cd(dataDir);

%%
manualMeta = struct();
manualMeta.nWells = 1;
meta = Metadata(dataDir, manualMeta);
meta.channelLabel = {'DAPI','CDX2','FOXA2','EPCAM'};

%filelist = dir(fullfile(dataDir,'*_FusionStitcher.ims')); % PICK FILES HERE
filelist = dir(fullfile(dataDir,'*_FusionStitcher.ims')); % PICK FILES HERE

%% go through list at fixed cross section positions

% just cross sections along x for now
index = 430;
RGBset = [1 2 4];
tol = [0.01 0.99; 0.01 0.99; 0.01 0.99];

for fi = 1%:numel(filelist)
    
    fname = filelist(fi).name;
    im = readStack(fname);
    
    % ------- z-correction: comment out for speed ---------
    % sub bg
    imbgsub = im;
    tic
    for zi = 1:size(im,4)
        for ci = 1:4
            slice = imgaussfilt(imbgsub(:,:,ci,zi), 2);
            bg = imopen(slice, strel('disk',50));
            imbgsub(:,:,ci,zi) = slice - bg;
        end
    end
    toc

    % correct intensity
    imcorr = imbgsub;
    f = [1.1 1.1 1.1, 1.1];
    tic
    for zi = 1:size(imcorr,4)
        for ci = 1:4
            imcorr(:,:,ci,zi) = imcorr(:,:,ci,zi)*(f(ci)^(zi-1));
        end
    end
    toc
    % ------- end comment out for speed ---------
    makeXsection(im, fname, meta, index, RGBset, [], tol);
end


%% alternative: read on file, then make multiple cross sections

fi = 1;
fname = filelist(fi).name;
im = readStack(fname);

%%
index = 1000;
RGBset = [2 3 4];
nucChannel  = 1;
 
figure,
makeXsection(imbgsub, fname, meta, index, RGBset, [], tol);


%% function definition:

function makeXsection(im, fname, meta, index, RGBset, nucChannel, tol)
    % im: image stack
    % meta: metadata
    % index : y index to make cross section at
    % RGBset : channels to use for RGB image

    lw = 10;

    Ilim = {};
    MIP = max(im,[],4);

    MIPadj = MIP; 
    for ci = 1:meta.nChannels
        Ilim{ci} = stretchlim(MIP(:,:,ci), tol(ci,:));
        MIPadj(:,:,ci) = imadjust(MIP(:,:,ci),Ilim{ci}); 
    end

    t = round(lw/2);
    MIPadj(index-t:index+t,:,:) = max(MIPadj(:));

    newnslices = round(meta.zres / meta.xres)*size(im,4);
    scaled = zeros([newnslices size(im,2)],'uint16');
    for ci = 1:meta.nChannels
        scaled(:,:,ci) = flipud(imresize(squeeze(im(index,:,ci,:)), [size(im,2) newnslices])');
        scaled(:,:,ci) = imadjust(scaled(:,:,ci),Ilim{ci});
    end

    %----- display

    figure,
    xyim = MIPadj(:,:,RGBset);
    if ~isempty(nucChannel)
        xyimnuc = MIPadj(:,:,nucChannel);
        for i = 1:3
           xyim(:,:,i) =  xyim(:,:,i) + xyimnuc;
        end
    end
    imshow(xyim,[])
    imwrite(xyim, [fname(1:end-4) '_' meta.channelLabel{RGBset} '_MIP.png']);

    figure,
    for ci = 1:meta.nChannels

        subplot_tight(meta.nChannels + 1, 1, ci)
        imshow(scaled(:,:,ci),[])
        title(meta.channelLabel{ci})
    end

    subplot_tight(meta.nChannels + 1, 1, meta.nChannels + 1)
    xzim = scaled(:,:,RGBset);
    if ~isempty(nucChannel)
        xzimnuc = scaled(:,:,nucChannel);
        for i = 1:3
           xzim(:,:,i) =  xzim(:,:,i) + xzimnuc;
        end
    end
    imshow(xzim,[])

    str = [];
    for i = 1:3
        str = [str meta.channelLabel{RGBset(i)} '; '];
    end
    title(str);
    saveas(gcf, [fname(1:end-4) '_xsect_' num2str(index) '.png']);
end