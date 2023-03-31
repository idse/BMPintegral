function [MIPtot, MIPidxtot] = makeMIP_Andor(inputdir, position, channel, outputdir, saveidx, tmax, type)

    % this is 180611 version, newest version has bug in concatenating 
    % separated time series because of added z-slice functionality 
    % choosing not to debug

    % turn off annoying BF warning
    warning('off','BF:lowJavaMemory');
    
    % read metadata
    meta = MetadataAndor(inputdir);
    filenameFormat = meta.filename;
    s = strsplit(meta.filename,'_f');
    if numel(s) < 2
        s = strsplit(meta.filename,'_p');
    end
    barefname = s{1};

    % some input checking
    if ~exist('saveidx','var')
        saveidx = false;
    end
    if ~exist('outputdir','var')
        outputdir = inputdir;
    end
    if ~exist(outputdir,'dir')
        mkdir(outputdir);
    end
    if isempty(strfind(filenameFormat,'_w'))
        warning('this will go twice as fast if you export channels separately');
    end
    ci = channel;
    if ci+1 > meta.nChannels 
        error('channel exceeds number of channels');
    end
    pi = position;
    if pi+1 > meta.nPositions
        error('position exceeds number of positions');
    end
    
    if ~exist('tmax','var')
        tmax = Inf;
    end
    
    % deal with the case where z is separated
    if ~isempty(strfind(filenameFormat,'_z'))
        fileZmax = meta.nZslices;
        zslices = {};
    else
        fileZmax = 1;
    end

    MIP = {};
    MIPidx = {};

    fileTmax = ceil(meta.nTime/meta.tPerFile) - 1;
    
    for ti = 0:fileTmax

        if ti*meta.tPerFile < tmax

            if (ti + 1)*meta.tPerFile > tmax
                tmaxinfile = tmax - ti*meta.tPerFile;
            else
                tmaxinfile = meta.tPerFile;
            end

            % read the data 
            % (if all time points are in one file there is no _t part)
            for zi = 0:fileZmax-1

                if ~isempty(strfind(filenameFormat,'_w'))
                    if fileTmax > 0 
                        fname = sprintf(filenameFormat, pi, ti, ci);
                    else
                        fname = sprintf(filenameFormat, pi, ci);
                    end
                else
                    % this is assuming either split z or t, not both
                    if fileTmax > 0 
                        fname = sprintf(filenameFormat, pi, ti);
                    elseif fileZmax > 1
                        fname = sprintf(filenameFormat, pi, zi);
                    else
                        fname = sprintf(filenameFormat, pi);
                    end
                end

                tic
                zslices{zi+1} = readStack(fullfile(inputdir,fname));
                toc
            end
            stack = cat(3, zslices{:});

            % make MIP
            if strcmp(type,'MIP')

                if ~isempty(strfind(filenameFormat,'_w'))
                    [MIP{ti+1},MIPidx{ti+1}] = max(stack(:,:,:,1,1:tmaxinfile),[],3);
                else
                    [MIP{ti+1},MIPidx{ti+1}] = max(stack(:,:,:,channel+1,1:tmaxinfile),[],3);
                end
                MIPidx{ti+1} = uint16(MIPidx{ti+1});

            % or SIP
            elseif strcmp(type,'SIP')

                if ~isempty(strfind(filenameFormat,'_w'))
                    MIP{ti+1} = uint16(sum(stack(:,:,:,1,1:tmaxinfile),3));
                else
                    MIP{ti+1} = uint16(sum(stack(:,:,:,channel+1,1:tmaxinfile),3));
                end
                warning('SIP cast to uint16 bc imwrite doesnt take uint32, should fix this'); 
            end
            clear stack; 
        end
    end
    
    % combine time series from different files
    MIPtot = cat(5,MIP{:});
    MIPidxtot = cat(5,MIPidx{:});
    clear MIP MIPidx;

    % save result
    %-------------

    % MIP/SIP
    fname = fullfile(outputdir, sprintf([barefname '_' type '_p%.4d_w%.4d.tif'],pi,ci));
    if exist(fname,'file')
        delete(fname);
    end
    % I HAD PROBLEMS WITH BFSAVE
    %bfsave(MIPtot,fname, 'dimensionOrder', 'XYZCT');
    MIPtot = squeeze(MIPtot);
    imwrite(MIPtot(:,:,1), fname);
    for i = 2:size(MIPtot,3)
        imwrite(MIPtot(:,:,i), fname,'WriteMode','Append');
    end

    % MIP index
    if saveidx
        fname = fullfile(outputdir, sprintf([barefname '_MIPidx_p%.4d_w%.4d.tif'],pi,ci));
        if exist(fname,'file')
            delete(fname);
        end
        %bfsave(MIPidxtot,fname, 'dimensionOrder', 'XYZCT');
        MIPidxtot = squeeze(MIPidxtot);
        imwrite(MIPidxtot(:,:,1), fname);
        for i = 2:size(MIPidxtot,3)
            imwrite(MIPidxtot(:,:,i), fname, 'WriteMode','Append');
        end
    end
    
    warning('on','BF:lowJavaMemory');
end