function upperleft = stitchedPreviews(dataDir, meta, options, MIPoptions)
% stitchedPreviewsFusion(dataDir, meta, options)
% options: struct with fields
% - margin
% - montages : either number of conditions (regular) or number of colonies
%               / grids in case of micropatterns / montages 
% - ss: subsampling (size reduction) factor for previews
% - saveStitchFullSize
% - FusionMontage : boolean if using the montage from Fusion
% - dichroic correction channels: 
%     which channels to correct for the CYR dichroic relative to BGRI (count from 0)

if ~isfield(options, 'saveStitchFullSize')
    options.saveStitchFullSize = false;
end
if ~isfield(options, 'FusionMontage')
    options.FusionMontage = false;
end
if ~isfield(options, 'dichroicCorrectionChannels')
    options.dichroicCorrectionChannels = [];
end
if ~isfield(options, 'type')
	options.type = 'MIP';
end
if ~isfield(options, 'margin')
    margin = 0;
else
    margin = options.margin;
end
if ~isfield(options,'montages')
    options.montages = 1:meta.nWells;
end
if ~isfield(options,'stitchStacks')
    options.stitchStacks = false;
end
if ~isfield(options,'videoOptions')
    videoOptions = struct('format','mp4', 'FrameRate', 5);
end
if ~isfield(options,'nucChannel')
    nucChannel = 0;
else
    nucChannel = options.nucChannel;
end
if ~isfield(options,'stitchMethod')
    options.stitchMethod = 'maximum';
end

% load dichroic correction 
dichroicCorr = load('dichroic_transform.mat');

% subsampling factor
if ~isfield(options, 'ss')
    if meta.xSize <= 1024
        ss = 2;
    else
        ss = 4;
    end
else
    ss = options.ss;
end
fullmargin = margin;
margin = round(margin/ss);

MIPfiles = dir(fullfile(dataDir,'MIP',['*' options.type '_*tif']));
s = strsplit(MIPfiles(1).name,['_' options.type]);
barefname = s{1};

pixelOverlap = round(meta.xSize*meta.montageOverlap/100);
posPerCondition = meta.posPerCondition;

if ~isempty(meta.montageGridSize)
    gridSize = meta.montageGridSize;
else % when not actually stitching
    gridSize = [round(posPerCondition/2) 2];
    meta.montageOverlap = 0;
end

% estimated size of stitched image
stitchedSize = round([  gridSize(1) - meta.montageOverlap*(gridSize(1)-1)/100,... 
                        gridSize(2) - meta.montageOverlap*(gridSize(2)-1)/100]*1024);
                    
%-------- colony names for Fusion ---------------------
if options.FusionMontage
    
    listing = dir(fullfile(dataDir,'MIP','*_MIP_p0000_w0000.tif'));
    exclude = contains({listing.name},'stitched'); 
    listing = listing(~exclude);
    
    colbarefnames = {};
    for i = 1:numel(listing)
        s = strsplit(listing(i).name,'_MIP_');
        colbarefnames{i} = s{1};
    end
    options.montages = 1:numel(colbarefnames);
    posPerCondition = prod(meta.montageGridSize);
    
    if posPerCondition == 9
        conditionPositionFusion = [7 2 9; 4 5 6; 1 8 3];
    elseif posPerCondition == 4
        conditionPositionFusion = [3 4; 1 2];
    else
        conditionPositionFusion = flipud(conditionPositionFusion');
    end
end
%------------------------------------------------------

upperleft = {};

fileID = fopen(fullfile(dataDir,'filenames.txt'),'w');
for montagenr = options.montages
    
    if options.FusionMontage
        colbarefname = colbarefnames{montagenr};
        conditionPositions = 1:posPerCondition;
    else
        colbarefname = barefname;
        conditionPositions = posPerCondition*(montagenr-1)+1:posPerCondition*montagenr;
    end
    fnameformat = [colbarefname '_' options.type '_p%.4d_w%.4d.tif'];
    disp(['file : ' colbarefname]);
    
    str1 = sprintf('stitched_MIP_p%.4d',montagenr-1);
    fprintf(fileID,(strcat(colbarefname, " -> ", str1, "\n")));
    
    channelorder = [nucChannel+1, 1:nucChannel,...
        nucChannel+2:numel(MIPoptions.channels)];
    for cii = channelorder
        
        
        ci = MIPoptions.channels(cii);
        
        disp(['-------------processing channel ' num2str(ci) '--------']);
        
        imgs = {};
        idximgs = {};
        rawimgs = {};
        if isfield(MIPoptions, 'tmax') && ~isempty(MIPoptions.tmax)
            tmax = MIPoptions.tmax;
        else
            tmax = meta.nTime;
        end
        
        for ti = 1:tmax
            
            disp(['processing time ' num2str(ti)]);

            %-----------------------------------------------------
            % read all the data
            %-----------------------------------------------------
            
            for pi = conditionPositions

                % gridSize 1 and 2 may be swapped, I have no way of knowing right now
                if options.FusionMontage
                    [i,j] = ind2sub(gridSize, conditionPositionFusion(pi));
                else
                    [i,j] = ind2sub(gridSize, pi - conditionPositions(1) + 1);
                end
                
                disp(['reading MIP ' num2str(pi)]);
                %----------------------------------------------------------------------
                fname = fullfile(dataDir,'MIP',sprintf(fnameformat,pi-1, ci));
                if exist(fname,'file') % to deal with missing files, just add black image
                    imgs{j,i} = imread(fname,ti);
                    if ismember(ci, options.dichroicCorrectionChannels)
                        disp('applying dichroic correction');
                        imgs{j,i} = imwarp(imgs{j,i}, dichroicCorr.tform, 'OutputView',imref2d(size(imgs{j,i})));
                    end
                else
                    warning(['channel ' num2str(ci) ' has MIP missing for position ' num2str(pi) ' (did you make MIP for all channels?)']);
                    imgs{j,i} = zeros([1024 1024],'uint16'); 
                end
                
                if MIPoptions.saveidx(cii) && options.stitchStacks
                    
                    disp(['reading MIPidx ' num2str(pi)]);
                    %------------------------------------------------------------------
                    mipidxFnameformat = [colbarefname '_MIPidx_p%.4d_w%.4d.tif'];
                    fname = fullfile(dataDir,'MIP',sprintf(mipidxFnameformat,pi-1, ci));
                    idximgs{j,i} = imread(fname,ti);
                end

                disp(['reading raw data for stitched stack panel ' num2str(pi)]);
                %----------------------------------------------------------------------
                if options.stitchStacks

                    % it seems Andor zero-pads depending on total
                    % this is a dumb solution up to 100 positions, do nicely later
                    rawfname = fullfile(dataDir, sprintf([colbarefname '_F%.2d.ims'],pi-1));
                    if ~exist(rawfname, 'file')
                        rawfname = fullfile(dataDir, sprintf([colbarefname '_F%.1d.ims'],pi-1));
                    end
                    if ~exist(rawfname, 'file')
                        rawfname = fullfile(dataDir, sprintf([colbarefname '_F%.3d.ims'],pi-1));
                    end
                    
                    if exist(rawfname,'file') % to deal with missing files, just add black image
                        
%                         stack = squeeze(readStack(rawfname, ci+1,ti));
                        stack = squeeze(readFrame(rawfname,ci+1,ti));
%                         stack = stack(:,:,:,end);
                        for zi = 1:meta.nZslices    
                            rawimgs{j,i,zi} = stack(:,:,zi);
                            if ismember(ci, options.dichroicCorrectionChannels)
                                disp('applying dichroic correction');
                                rawimgs{j,i,zi} = imwarp(rawimgs{j,i,zi}, dichroicCorr.tform, 'OutputView',imref2d(size(rawimgs{j,i,zi})));
                            end
                        end
                    else
                        error(['channel ' num2str(ci) ' has raw data missing for position ' num2str(pi)]);
                    end
                end
            end
            
            %-----------------------------------------------------
            % stitch together 
            %-----------------------------------------------------

            % stitch only the first channel
            if cii == nucChannel + 1
                
                % detect size of overlap in pixels
                dim = 2; 
                maxPixelOverlap = round(size(imgs{1,1}, dim)/4);
                offset = cross_corr(imgs{1,1},imgs{1,2}, dim, maxPixelOverlap);

                disp(['stitching based on channel ' num2str(ci)]);
                if ~isempty(pixelOverlap) %% && ti == 1
                    
                    if abs(offset/pixelOverlap - 1) > 0.25
                        warning('detected pixel overlap very different from nominal value, sticking with nominal value');
                    else
                        pixelOverlap = offset;
                    end
                    % get register positions of upper left corner
                    upperleft{montagenr}{ti} = registerImageGrid_v3(imgs, pixelOverlap);
                else %if ti == 1 && isempty(pixelOverlap)
                    for pi = conditionPositions
                        
                        if options.FusionMontage
                            [i,j] = ind2sgrub(gridSize, conditionPositionFusion(pi));
                        else
                            % don't understand anymore what is going on with index order but this seemed to work
                            % 200417 changed it again from [j,i] to [i,j],
                            % something must be flipped in a context
                            % dependent way, figure out later
                            [i,j] = ind2sub(gridSize, pi - conditionPositions(1) + 1);
                        end
                        upperleft{montagenr}{ti}{j,i} = [1+(j-1)*(meta.ySize + 50), 1+(i-1)*(meta.xSize + 50)];
                    end
                end
            end
            
            % CAUTION : upperleft = [y x]
            if strcmp(options.stitchMethod,'maximum')
                [stitched, upperleft{montagenr}{ti}] =...
                    stitchImageGrid(upperleft{montagenr}{ti}, imgs);
            elseif strcmp(options.stitchMethod,'weightedAvg')
                [stitched, upperleft{montagenr}{ti}] =...
                    stitchImageGridDistWeighted_v2(upperleft{montagenr}{ti}, imgs);
            end
            
            if ~isempty(pixelOverlap)
                % trim to projected size
                tmp = zeros(stitchedSize,'uint16');
                yrange = 1:min(stitchedSize(1),size(stitched,1));
                xrange = 1:min(stitchedSize(2),size(stitched,2));
                tmp(yrange, xrange) = stitched(yrange, xrange);
                stitched = tmp;
            end
            
            % make clean preview (not for quantitative analysis)
            small = imfilter(stitched,ones(ss)/ss^2);
            small = small(1:ss:end,1:ss:end);
            if ti == 1
                preview = zeros([size(small) tmax],'uint16');
            end
            preview(1:size(small,1), 1:size(small,2), ti) = small;
            
            %-----------------------------------------------------
            % save fullsize stitched MIP
            %-----------------------------------------------------
            if options.saveStitchFullSize
                
                if ti==1 
                    fullsize = size(stitched);
                end

                fname = sprintf('stitched_MIP_p%.4d_w%.4d_t%.4d.tif',montagenr-1,ci,ti-1);
               
                tmp = zeros(fullsize,'uint16');
                idx = min(fullsize, size(stitched));
                tmp(1:idx(1),1:idx(2)) = stitched(1:idx(1),1:idx(2));
                tmp = tmp(1+fullmargin:end-fullmargin, 1+fullmargin:end-fullmargin);
                
                imwrite(tmp, fullfile(dataDir,'MIP', fname));
            end
            
            %-----------------------------------------------------
            % save fullsize stitched stack
            %-----------------------------------------------------
            
            if options.stitchStacks
                
                stackfname = fullfile(dataDir,sprintf('stitched_p%.4d_w%.4d_t%.4d.tif', montagenr-1, ci, ti-1));
                
                disp(['writing stitched stack ' stackfname])
                for zi = 1:meta.nZslices
                    
                    if strcmp(options.stitchMethod,'maximum')
                        [stitchedSlice, upperleft{montagenr}{ti}] =...
                            stitchImageGrid(upperleft{montagenr}{ti}, rawimgs(:,:,zi));
                    elseif strcmp(options.stitchMethod,'weightedAvg')
                        [stitchedSlice, upperleft{montagenr}{ti}] =...
                            stitchImageGridWeightedAvg(upperleft{montagenr}{ti}, rawimgs(:,:,zi));
                    end

%                     [stitchedSlice, upperleft{montagenr}{ti}] = stitchImageGrid(upperleft{montagenr}{ti}, rawimgs(:,:,zi));
                    
                    % trim to projected size
                    tmp = zeros(stitchedSize,'uint16');
                    yrange = 1:min(stitchedSize(1),size(stitchedSlice,1));
                    xrange = 1:min(stitchedSize(2),size(stitchedSlice,2));
                    tmp(yrange, xrange) = stitchedSlice(yrange, xrange);
                    stitchedSlice = tmp;
                    
                    if zi == 1
                        imwrite(stitchedSlice, stackfname);
                    else
                        imwrite(stitchedSlice, stackfname, 'WriteMode','append');
                    end
                end
                
                % MIPidx stitching and saving
                if MIPoptions.saveidx(cii)
                    if strcmp(options.stitchMethod,'maximum')
                        [stitchedidx, upperleft{montagenr}{ti}] =...
                            stitchImageGrid(upperleft{montagenr}{ti}, idximgs);
                    elseif strcmp(options.stitchMethod,'weightedAvg')
                        [stitchedidx, upperleft{montagenr}{ti}] =...
                            stitchImageGridWeightedAvg(upperleft{montagenr}{ti}, idximgs);
                        stitchedidx(stitchedidx > meta.nZslices) = meta.nZslices;
                        stitchedidx(stitchedidx == 0) = 1;
                    end
%                     [stitchedidx, upperleft{montagenr}{ti}] = stitchImageGrid(upperleft{montagenr}{ti}, idximgs);
                    tmp(yrange, xrange) = stitchedidx(yrange, xrange);
                    stitchedidx(yrange, xrange) = stitchedidx(yrange, xrange);
                    fname = sprintf('stitched_MIPidx_p%.4d_w%.4d_t%.4d.tif',montagenr-1,ci,ti-1);
                    imwrite(tmp, fullfile(dataDir,'MIP', fname));
                end
            end
        end
        
        %-----------------------------------------------------
        % save preview avi + initial and final frame
        %-----------------------------------------------------
            
        % set lookup table
        MIPinTime = max(preview,[],3);              
        Ilim = stitchedlim(MIPinTime);

        % scale lim back to 16 bit range
        Imin = double(min(MIPinTime(MIPinTime>0))); 
        Imax = round(Ilim(2)*(2^16-1));
        preview = mat2gray(preview, [Imin Imax]);
        preview = uint16((2^16-1)*preview);

        % crop
        preview = preview(1+margin:end-margin, 1+margin:end-margin,:); 

        bareprevfname = sprintf('stichedPreview_p%.4d_w%.4d',montagenr-1, ci);
        imwrite(mat2gray(preview(:,:,1)),fullfile(dataDir, [bareprevfname '_initial.jpg']));
        if tmax > 1

            imwrite(mat2gray(preview(:,:,tmax)),fullfile(dataDir, [bareprevfname '_final.jpg']));

            % save as compressed video
            if strcmp(videoOptions.format,'avi')
                v = VideoWriter(fullfile(dataDir,[bareprevfname '.avi']),'Uncompressed AVI');
            elseif strcmp(videoOptions.format,'mp4')
                v = VideoWriter(fullfile(dataDir,[bareprevfname '.mp4']),'MPEG-4');
            else
                error('unknown video format');
            end
            v.FrameRate = videoOptions.FrameRate;
            v.Quality = 100;
            open(v)
            for ti = 1:tmax
                writeVideo(v,mat2gray(preview(:,:,ti)))
            end
            close(v);
        end
    end
end
fclose(fileID);
save(fullfile(dataDir,'upperleft'), 'upperleft');
end