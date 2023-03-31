function batchMIP(inputdir, outputdir, opts)

    % MIPoptions: struct with fields
    % channels, saveidx, tmax, positions

    % read metadata
    meta = Metadata(inputdir);
    
    % get file list
    extension = meta.filenameFormat(end-3:end);
    files = dir(fullfile(inputdir,['*' extension]));
    
    if ~isfield(opts,'tmax') || isempty(opts.tmax)
        opts.tmax = meta.nTime;
    end
    if ~isfield(opts,'positions')
        opts.positions = [];
    end
    if ~isfield(opts,'zrange')
        opts.zrange = [];
    end

    for i = 1:numel(files)
        
        filename = files(i).name;
        [splitfname, matches] = strsplit(filename,{'_F[0-9][0-9]*', ['\' extension]},...% F->. for any letter
                                   'DelimiterType','RegularExpression',...
                                   'CollapseDelimiters',false);
        barefname = splitfname{1};
        pi = str2num(matches{1}(3:end));
        
        if any(opts.positions==pi) || isempty(opts.positions)
        
            for ci = 1:numel(opts.channels)
                disp(['processing position ' num2str(pi) ', channel ' num2str(opts.channels(ci))]);
                makeMIP_bioformats(fullfile(inputdir,filename), barefname,...
                    pi, opts.channels(ci), outputdir, opts.saveidx(ci), opts.tmax, opts.zrange);
            end
        end
    end
end