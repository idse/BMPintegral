function batchMIP_Andor(inputdir, outputdir, channels, saveidx, tmax, positions, type)

    meta = MetadataAndor(inputdir);

    if ~exist('positions','var')
        positions = 1:meta.nPositions;
    end

    if ~exist('tmax','var')
        tmax = meta.nTime;
    end
    
    if ~exist('type','var')
        type = 'MIP';
    end
    
    for ci = 1:numel(channels)
        for pi = positions
            disp(['processing position ' num2str(pi) ', channel ' num2str(channels(ci))]);
            makeMIP_Andor(inputdir, pi-1, channels(ci), outputdir, saveidx(ci), tmax, type);
        end
    end
end