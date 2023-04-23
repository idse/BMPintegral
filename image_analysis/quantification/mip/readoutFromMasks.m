function cellData = readoutFromMasks(imgs,masks,opts)

if nargin < 3
    opts = struct();
end

mn = size(imgs{1},[1,2]); m = mn(1); n = mn(2);

if ~isfield(opts,'nucChannel')
    nucChannel = 0;
else
    nucChannel = opts.nucChannel;
end
if ~isfield(opts,'cytoplasmicLevels')
    opts.cytoplasmicLevels = false;
end
if ~isfield(opts,'zmethod')
    opts.zmethod = 'median';
end
if ~isfield(opts,'readoutmethod')
    opts.readoutmethod = 'mean';
end
if ~isfield(opts,'useMIP')
    opts.useMIP = false;
end

nchannels = length(imgs);
ncells = length(masks.nucmask);
nz = size(imgs{1},3);

if opts.useMIP
    for ii = 1:length(imgs)
        imgs{ii} = max(imgs{ii},[],3);
    end
    if opts.cytoplasmicLevels
        cytmasks = cell(1,ncells);
        for ii = 1:ncells
            cytmasks{ii} = unique(cell2mat(masks.cytmask(:,1)));
        end
    end
end

%optionally pass segmentation of junk in one or more channels
%incorporate this into masks in the future?
if isfield(opts,'junkmask') && ~isempty(opts.junkmask)
    %for each channel, replace image values with NaN inside the junk mask,
    %and ignore NaN values when extracting image data later
    for ii = 1:length(imgs)
        if ~isempty(opts.junkmask{ii})
            jm = opts.junkmask{ii};
            if size(jm,3) == 1 && nz > 1
                %if mip was segmented for junk, apply the same mask to all
                %slices, otherwise could segment junk across z and keep
                %that info
                jm = repmat(jm,1,1,nz);
            end
            imgs{ii}(jm==1) = NaN;
        end
    end
end

nucim = imgs{nucChannel+1};
[~,MIPidx] = max(nucim,[],3); %calculate MIPidx
MIPidx(MIPidx == 0) = 1;

cellData = struct;
cellData.cytLevel = zeros([ncells nchannels]);
cellData.NCratio = zeros([ncells nchannels]);
cellData.nucLevel = zeros([ncells nchannels]);
cellData.nucZ = zeros([ncells 1]);
cellData.background = zeros([1 nchannels]);

%nuclear geometry
cellData.XY = masks.XY;
cellData.nucArea = cellfun(@numel,masks.nucmask)';
if isfield(masks.nucstats,'Orientation')
    cellData.nucOrientation = cat(1,masks.nucstats.Orientation);
    cellData.nucMajorAxis = cat(1,masks.nucstats.MajorAxisLength);
    cellData.nucMinorAxis = cat(1,masks.nucstats.MinorAxisLength);
    cellData.nucCircularity = cat(1,masks.nucstats.Circularity);
elseif isfield(masks.nucstats,'nucOrientation')
    cellData.nucOrientation = cat(1,masks.nucstats.nucOrientation);
    cellData.nucMajorAxis = cat(1,masks.nucstats.nucMajorAxis);
    cellData.nucMinorAxis = cat(1,masks.nucstats.nucMinorAxis);
    cellData.nucCircularity = cat(1,masks.nucstats.nucCircularity);
end

nucmask = false(m,n);
nucmask(cell2mat(masks.nucmask(:))) = true;
bgmask = masks.bgmask;

% for background subtraction, median/mode z-plane
if strcmp(opts.zmethod,'mode')
    zmed = mode(MIPidx(nucmask));
elseif strcmp(opts.zmethod,'median')
    zmed = round(median(MIPidx(nucmask)));
end

% read out nuclear and cytoplasmic levels
%-----------------------------------------
if ncells > 0
    for cii = 1:nchannels
        imc = imgs{cii};
        
        % current low-tech background subtraction:
        %-------------------------------------------------
        % mean value of segmented empty space in the image
        % or otherwise just min of image
        if ~isempty(bgmask)
            % if the background area is too small to be
            % reliable (< ~1% of total area), then leave as NaN

            imcZmed = imc(:,:,zmed);
            %dont use black buffer regions around the image for
            %background estimation
            bgmask = bgmask & (imcZmed > 0);
            if sum(bgmask,'all') > numel(bgmask)/100
                %use weighted average over z slices for background estimation
                vals = MIPidx(nucmask);
                indices = unique(vals);
                npoints = numel(vals);
                bg_est = 0;
                for jj = 1:length(indices)
                    ind = indices(jj);
                    imcZmed = imc(:,:,ind);
                    weight = sum(vals == ind)/npoints;
                    bg_est = bg_est + weight*mean(imcZmed(bgmask),'omitnan');
                end
                cellData.background(cii) = bg_est;
            else
                %not a good solution -> better to use previous info
                %could handle that but it is easier to do outside of this
                %function
                cellData.background(cii) = NaN;
            end
        else
            imcZmed = imc(:,:,zmed);
            cellData.background(cii) = min(imcZmed(imcZmed>0),[],'all');
        end
        
        for cellidx = 1:ncells
            
            nucPixIdx = masks.nucmask{cellidx};
            if strcmp(opts.zmethod,'median')
                zi = round(median(MIPidx(nucPixIdx)));
            elseif strcmp(opts.zmethod,'mode')
                zi = mode(MIPidx(nucPixIdx));
            end
            
            nucPixIdx = double(zi-1)*m*n + nucPixIdx;
            if strcmp(opts.readoutmethod,'mean')
                cellData.nucLevel(cellidx, cii) = mean(imc(nucPixIdx),'omitnan');
            elseif strcmp(opts.readoutmethod,'median')
                cellData.nucLevel(cellidx, cii) = median(imc(nucPixIdx),'omitnan');
            end
            cellData.nucZ(cellidx) = zi;

            if opts.cytoplasmicLevels
                if size(masks.cytmask,1) == nz
                    cytPixIdx = masks.cytmask{zi,cellidx};
                elseif size(masks.cytmask,1) == 1
                    cytPixIdx = masks.cytmask{cellidx};
                else
                    error('cytmask z size does not match number of z slices')
                end
                cytPixIdx = double(zi-1)*m*n + cytPixIdx;
                
                if strcmp(opts.readoutmethod,'mean')
                    cellData.cytLevel(cellidx, cii) = mean(imc(cytPixIdx),'omitnan');
                elseif strcmp(opts.readoutmethod,'median')
                    cellData.cytLevel(cellidx, cii) = median(imc(cytPixIdx),'omitnan');
                end
                
                C = cellData.cytLevel(cellidx, cii);
                N = cellData.nucLevel(cellidx, cii);
                bg = cellData.background(cii);
                cellData.NCratio(cellidx, cii) = (N-bg)./(C-bg);
            end
        end
    end
else
    warning('------------ NO CELLS ------------');
end



end