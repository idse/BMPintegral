function newCellData = extractDataMultiZ(dataDir, fname, meta, opts)

    %options for consolidation of cell info between z stacks
    if ~isfield(opts,'IoU')
        opts.IoU = 0.75;
    end
    if ~isfield(opts,'maxZsize')
        opts.maxZsize = 15; %um
    end
    if ~isfield(opts,'minZsize')
        opts.minZsize = 1; %um
    end
    if ~isfield(opts,'maxCentroidDistance')
        opts.maxCentroidDistance = 5; %um
    end
    %z method options:
        %read cell data in all channels based on the channel in which the
        %nuclear marker (usually DAPI) is brightest
        %read out the nuclear intensity in each channel in a given nucleus
        %in the z slice in which the intensity in that channel is highest
        %read out the intensity a fixed number of z slices away from the
        %slice in which the DAPI signal is brightest (account for chromatic
        %aberration)
    if ~isfield(opts,'Zmethod')
        opts.Zmethod = 'useNuclearMarker';
    end
    if ~isfield(opts,'LM_as_cell')
        opts.LM_as_cell = false;
    end
    
    if contains(fname,'stitched_p')
        prefix = fname(1:14);
    elseif contains(fname,'FusionStitcher')
        s = strsplit(fname,{'_FusionStitcher','.ims'});
        prefix = [s{:}];
    elseif contains(fname,'_Stitched') && contains(fname,'.tif')
        s = strsplit(fname,{'_Stitched','.tif'});
        prefix = [s{:}];
    elseif contains(fname,'.tif')
        s = strsplit(fname,{'.tif'});
        prefix = [s{:}];
    elseif contains(fname,'.nd2')
        s = strsplit(fname,{'.nd2'});
        prefix = [s{:}];
    else
        s = strsplit(fname,'.ims');
        prefix = [s{:}];
    end
    subDir = fullfile(dataDir,[prefix '_zslices']);
    save(fullfile(subDir,'meta.mat'),'meta')
    
    pattern = [prefix '_w%.4d_z%.4d'];
    
    cellData(meta.nZslices) = struct;
    ncells = zeros(meta.nZslices,1);
    
    %initialize bioformats reader
    r = bfGetReader(fullfile(dataDir, fname));
    r.setSeries(0);
    %compute MIPs during processing step
    MIPs = cell(1,meta.nChannels);
    
    for zi = 1:meta.nZslices
        
        fprintf('.')
        maskname = [fullfile(subDir,sprintf(pattern,meta.nucChannel,zi-1)),'_FinalSegmentation.png'];
        %read in the label matrix 
        mask = imread(maskname);
        imsize = size(mask);
        %binarize and remove the boundaries of objects so they can be
        %identified as connected components
        if ~islogical(mask) %if the mask is a label matrix
            nucmask = mask > 0 & ~boundarymask(mask);
        else %if the mask is already binarized
            nucmask = mask;
        end

        %make cytoplasmic masks
        nucmaskmarg = imdilate(nucmask,strel('disk',opts.cytoMargin));
        dilated = imdilate(nucmask, strel('disk',opts.cytoSize + opts.cytoMargin));
        basin = bwdist(dilated);
        basin = imimposemin(basin, nucmask);
        L = watershed(basin);

        nucstats = regionprops(L, 'PixelIdxList');
        cytCC = struct('PixelIdxList', {cat(1,{nucstats.PixelIdxList})});

        for cci = 1:numel(cytCC.PixelIdxList)
            CCPIL = cytCC.PixelIdxList{cci};
            CCPIL = CCPIL(dilated(CCPIL));    % exclude outside dilated nuclei
            CCPIL = CCPIL(~nucmaskmarg(CCPIL)); % exclude margin around nuclei
            cytCC.PixelIdxList{cci} = CCPIL;
        end

        % make it a proper CC struct
        cytCC.Connectivity = 8;
        cytCC.NumObjects = size(cytCC.PixelIdxList,2);
        cytCC.ImageSize = size(nucmask);
        cytstats = regionprops(cytCC,'Area');
%         cytmask = false(size(nucmask));
%         cytmask(cat(1,cytCC.PixelIdxList{:}))=true;

        cellData(zi).cytMaskArea = [cytstats.Area]';

        cellData(zi).cytLevel = zeros([numel(cytCC.PixelIdxList) meta.nChannels]);
        cellData(zi).NCratio = zeros([numel(cytCC.PixelIdxList) meta.nChannels]);

        nucstats = regionprops(nucmask, 'Area', 'Centroid','PixelIdxList',...
                                        'Orientation', 'MajorAxisLength','MinorAxisLength',...
                                        'Circularity');

        nCells = numel(nucstats);

        ncells(zi) = nCells;
        cellData(zi).XY = cat(1,nucstats.Centroid);

        % nuclear geometry
        cellData(zi).nucArea = cat(1,nucstats.Area);
        cellData(zi).nucOrientation = cat(1,nucstats.Orientation);
        cellData(zi).nucMajorAxis = cat(1,nucstats.MajorAxisLength);
        cellData(zi).nucMinorAxis = cat(1,nucstats.MinorAxisLength);
        cellData(zi).nucCircularity = cat(1,nucstats.Circularity);

        cellData(zi).PixelIdxList = {nucstats.PixelIdxList};

        cellData(zi).nucLevel = zeros([nCells meta.nChannels]);
        cellData(zi).nucZ = zeros([nCells 1]);

        cellData(zi).background = zeros([1 meta.nChannels]);

        for cii = 1:meta.nChannels
            %read image
            img = bfGetPlane(r, r.getIndex(zi-1,cii-1,0)+1);
            
            %MIP
            if isempty(MIPs{cii})
                MIPs{cii} = img;
            else
                MIPs{cii} = max(MIPs{cii},img);
            end

            L = 200;
%             bg = min([mean(img(1:L,1:L),'all'), mean(img(1:L,end-L:end),'all'),...
%                  mean(img(end-L:end,1:L),'all'), mean(img(end-L:end,end-L:end),'all')]);
             
             bgBlock = {img(1:L,1:L),img(1:L,end-L:end),...
                 img(end-L:end,1:L),img(end-L:end,end-L:end)};
             bg = min(cellfun(@(x) mean(x(x>0),'all'),bgBlock));

            %bg = mean(img(yl,xl),'all');
            
            cellData(zi).background(cii) = bg;
            %get image intensity information for each cell in this channel and
            %slice
            for cellidx = 1:nCells
                %nuclear level
                nucPixIdx = nucstats(cellidx).PixelIdxList;
                cellData(zi).nucLevel(cellidx, cii) = mean(img(nucPixIdx));
                cellData(zi).nucZ(cellidx) = zi;

                %cytoplasmic level and NC ratio
                cytPixIdx = cytCC.PixelIdxList{cellidx};
                cellData(zi).cytLevel(cellidx, cii) = mean(img(cytPixIdx));
                C = cellData(zi).cytLevel(cellidx, cii);
                N = cellData(zi).nucLevel(cellidx, cii);
                bg = cellData(zi).background(cii);
                cellData(zi).NCratio(cellidx, cii) = (N-bg)./(C-bg);
            end
        end
    end
    fprintf('\n')
    
    [newCellData, chain] = combineZsegmentation(cellData, meta, opts);    

    %make 3D label matrix
    LM = cell(1,1,meta.nZslices);
    for zi = 1:meta.nZslices
        LM{zi} = zeros(imsize);
    end
    
    for ii = 1:length(chain)
        for zi = 1:size(chain{ii},1)
            z = chain{ii}(zi,1); %z index
            c = chain{ii}(zi,2); %cell index
            %pixel indices of mask for cell c in stack z
            PixelIdxList = cellData(z).PixelIdxList{c};
            LM{z}(PixelIdxList) = ii;
        end
    end
    
    LM = cell2mat(LM);
    
    %for some datasets, LM is very large (generally when there are >20 z
    %slices), and we may want to save memory by converting to a sparse
    %array; 3d sparse arrays aren't handled, so use a cell of 2D sparse
    %slices
    %this also gets around matlab's limit on variable size in .mat files,
    %but we could alternatively do that by specifying -v7.3 in the save function
    if opts.LM_as_cell
        LMcell = cell(1,1,size(LM,3));
        for zi = 1:size(LM,3)
            LMcell{zi} = sparse(LM(:,:,zi));
        end
        LM = LMcell;
    end
    
    savename = fullfile(dataDir,[fname(1:end-4),'_LabelMatrix3D.mat']);
    save(savename,'LM')
    
    
    %visualize it all together on MIP
    MIP = MIPs{meta.nucChannel+1};
    MIP = im2double(imadjust(MIP,stitchedlim(MIP)));
    
    newseg = struct;
    
    for jj = 1:length(chain)
        cellinfo = chain{jj};
        [~, I] = max(cellinfo(:,3));
        zi = cellinfo(I,1);
        ci = cellinfo(I,2);
        newseg(jj).z = zi;
        newseg(jj).PixelIdxList = cellData(zi).PixelIdxList{ci};
    end
    
    ncomps = length(newseg);
    zs = [newseg.z]';
    colors = zeros(ncomps,3);
    for ci = 1:9000:ncomps
        idxs = ci:min(ncomps,ci+8999);
        colors(idxs,:) = distinguishable_colors(length(idxs),'k');
    end
    
    im1 = zeros(size(img)); im2 = im1; im3 = im1;
    for zi = meta.nZslices:-1:1
        idxs = find(zs == zi);
        idxlists = {newseg(idxs).PixelIdxList};
        for li = 1:length(idxs)
            im1(idxlists{li}) = colors(idxs(li),1);
            im2(idxlists{li}) = colors(idxs(li),2);
            im3(idxlists{li}) = colors(idxs(li),3);
        end
    end
    
    overlay = 0.5*cat(3,im1,im2,im3) + 0.5*repmat(MIP,1,1,3);
    savename = fullfile(dataDir,[fname(1:end-4),'_SegOverlay.jpeg']);
    imwrite(overlay,savename)
    imshow(overlay)
        
end