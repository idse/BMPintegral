function [img, omeMeta] = readFrame(fullfname, channels, time) 
    % read the data from a single multichannel stack
    % same as readStack but only reads a single time point instead of
    % reading in the full stack from 1 to tmax
            
    % load the Bio-Formats library into the MATLAB environment
    autoloadBioFormats = 1;
    status = bfCheckJavaPath(autoloadBioFormats);
    assert(status, ['Missing Bio-Formats library. Either add loci_tools.jar '...
        'to the static Java path or add it to the Matlab path.']);

    % initialize logging
    loci.common.DebugTools.enableLogging('INFO');
    
    % create bioformats reader for file
    disp(fullfname);
    
    tic
    r = bfGetReader(fullfname);
    if ~exist('channels','var') || isempty(channels)
        channels = 1:r.getSizeC();
    end
    if ~exist('time','var')
        time = 1;%r.getSizeT();
    end
    img = zeros([r.getSizeY(), r.getSizeX(), numel(channels), r.getSizeZ()], 'uint16');
    for cii = 1:numel(channels)
        for zi = 1:r.getSizeZ()
            fprintf('.');
            readidx = r.getIndex(zi-1,channels(cii)-1,time-1)+1;
            disp(readidx)
            img(:,:,cii,zi) = bfGetPlane(r, readidx);
        end
    end
    omeMeta = r.getMetadataStore();
    r.close();
    toc
    %img = squeeze(img); % DON'T DO THIS, makeMIP needs xyczt even if c is
    %singleton because of the channels parameter
end