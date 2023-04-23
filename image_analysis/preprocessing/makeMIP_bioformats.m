function MIP = makeMIP_bioformats(...
    filename, barefname, position, channel, outputdir, saveidx, tmax, zrange)

    r = bfGetReader(filename);
    dataDir = fileparts(filename);

    % some input checking
    if ~exist('saveidx','var')
        saveidx = false;
    end
    if ~exist('outputdir','var') || isempty(outputdir)
        outputdir = fullfile(dataDir,'MIP');
    end
    if ~exist(outputdir,'dir')
        mkdir(outputdir);
    end
    
    % we only want to loop over series for lif files
    % because vsi files have the same picture at different resolution as
    % the series
    if strcmp(filename(end-3:end),'.lif') || strcmp(filename(end-3:end),'.nd2')
        N = r.getSeriesCount;
    else
        N = 1;
    end
    
    for si = 1:N
        
        disp(['reading series ' num2str(si)]);
        r.setSeries(si-1);
        
        ci = channel;
        if ci+1 > r.getSizeC()
            error('channel exceeds number of channels');
        end
        if ~exist('tmax','var') || isempty(tmax)
            tmax = r.getSizeT();
        else
            tmax = min(r.getSizeT(), tmax);
        end
        if ~exist('zrange','var') || isempty(zrange)
            zrange = 1:r.getSizeZ();
        end

        MIP = zeros([r.getSizeY() r.getSizeX() tmax], 'uint16');
        MIPidx = zeros([r.getSizeY() r.getSizeX() tmax], 'uint8');

        for ti = 1:tmax
            if numel(zrange) > 1
                im = zeros([r.getSizeY() r.getSizeX() numel(zrange)]);
                for zi = 1:numel(zrange)
                    im(:,:,zi) = bfGetPlane(r, r.getIndex(zrange(zi)-1,ci,ti-1)+1);
                end
                [MIP(:,:,ti), MIPidx(:,:,ti)] = max(im,[],3);
            else
                MIP(:,:,ti) = bfGetPlane(r, r.getIndex(0,ci,ti-1)+1);
            end
            fprintf('.');
            if mod(ti,60)==0
                fprintf('\n');
            end
        end
        fprintf('\n');

        % save result
        %-------------

        % MIP
        
        if r.getSeriesCount > 1 && (strcmp(filename(end-3:end),'.lif') || strcmp(filename(end-3:end),'.nd2'))
            pi = si-1;
            warning('multiple series lif, assuming all positions are in single file');
        else
            pi = position;
        end

        if ~isempty(pi)
            fname = fullfile(outputdir, sprintf([barefname '_MIP_p%.4d_w%.4d.tif'],pi,ci));
            idxfname = fullfile(outputdir, sprintf([barefname '_MIPidx_p%.4d_w%.4d.tif'],pi,ci));
        else
            fname = fullfile(outputdir, sprintf([barefname '_MIP_w%.4d.tif'],ci));
            idxfname = fullfile(outputdir, sprintf([barefname '_MIPidx_w%.4d.tif'],ci));
        end
        if exist(fname,'file')
            delete(fname);
        end

        MIP = squeeze(MIP);
        disp('saving MIP');
        imwrite(MIP(:,:,1), fname);
        fprintf('.');
        for i = 2:size(MIP,3)
            fprintf('.');
            if mod(i,60)==0
                fprintf('\n');
            end
            imwrite(MIP(:,:,i), fname,'WriteMode','Append');
        end
        fprintf('\n');

        if saveidx
            disp('saving MIPidx');
            MIPidx = squeeze(MIPidx);
            fprintf('.');
            imwrite(MIPidx(:,:,1), idxfname);
            for i = 2:size(MIPidx,3)
                fprintf('.');
                if mod(i,60)==0
                    fprintf('\n');
                end
                imwrite(MIPidx(:,:,i), idxfname,'WriteMode','Append');
            end
            fprintf('\n');
        end
    end
	r.close();
end