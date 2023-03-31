function newNuclearMask = nuclearCleanup(nuclearMask, options)
    % clean up nuclear mask
    %
    % newNuclearMask = nuclearCleanup(nuclearMask, options)
    %
    % nuclearMask:      input binary mask of nuclei
    %
    % newNucleiMask:    clean mask
    %    
    % options:          structure with fields:
    %
    % -minArea          remove junk smaller than this (in pixels)
    % -openSize         disk radius for imopen strel
    %
    % -separateFused    boolean
    % -clearBorder      boolean
    % -clearMinSolidity delete objects that have solidity smaller than this
    % value
    %
    % options for separateFusedNuclei:
    %
    % -minAreaStd:      only objects with A > mean(A) + minAreaStd*std(A)
    %                   can be considered fused (default 1)
    % -minFusedArea     alternative to minAreaStd that is more robus to large numbers of
    %                   fusions but requires more thought
    % -minSolidity:     only objects with solidity less than this can be
    %                   considered fused (default 0.95)
    %                   NOTE: this part is computationally expensive
    %                   set value <= 0 to turn off and speed up
    % -erodeSize        in units of mean radius, default 1
    
    
    % ---------------------
    % Idse Heemskerk, 2016
    % ---------------------

    if ~exist('options','var')
        options = struct();
    end
    if ~isfield(options,'minArea')
        options.minArea = 100;
    end
    if ~isfield(options,'openSize')
        options.openSize = 5;
    end
    if ~isfield(options,'separateFused')
        options.separateFused = true;
    end
    if ~isfield(options,'clearMinSolidity')
        options.clearMinSolidity = 0;
    end
    if ~isfield(options,'fillholes')
        options.fillholes = true;
    end
    
    CC = bwconncomp(nuclearMask);
    if options.clearMinSolidity > 0
        stats = regionprops(CC, 'ConvexArea', 'Area');
        convexArea = [stats.ConvexArea];
    else
        stats = regionprops(CC, 'Area');
    end
    area = [stats.Area];

    if options.clearMinSolidity > 0
        deletables = area./convexArea < options.clearMinSolidity;
        sublist = CC.PixelIdxList(deletables);
        nuclearMask(cat(1,sublist{:})) = false;
    end
    
    nuclearMask = bwareaopen(nuclearMask, options.minArea/5);
    nuclearMask = imopen(nuclearMask, strel('disk',options.openSize));
    % fill smaller holes that can appear in nuclear segmentation:
    nuclearMask = ~bwareaopen(~nuclearMask,options.minArea/5); 
    
    if options.separateFused && sum(nuclearMask(:))>0
        [nuclearMask, fusedMask] = separateFusedNuclei(nuclearMask,options);
    end
    if options.clearBorder 
        nuclearMask = imclearborder(nuclearMask);
    end
%    figure, imshow(fusedMask)
    if options.fillholes
        nuclearMask = imfill(nuclearMask,'holes');
    end
    nuclearMask = bwareaopen(nuclearMask, options.minArea);

    newNuclearMask = nuclearMask;
end