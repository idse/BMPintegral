classdef MetadataAndor < Metadata
    % metadata with additional properties for Andor IQ output
    
    % ---------------------
    % Idse Heemskerk, 2016
    % ---------------------
    
    %   Andor only file info:
    
    properties

        tPerFile            % number of time points per file
    end
    
    methods
        
        function this = MetadataAndor(dataDir)

            this = this@Metadata();
            
            if nargin == 1
                this = this.read(dataDir);
            end

            % always true for the Andor:
            this.xSize = 1024;
            this.ySize = 1024;
        end
        
        function this = read(this, dataDir)
            
            rawMeta = struct();
            
            % get metadata filename
            listing = dir(fullfile(dataDir,'*.txt'));
            if isempty(listing)
               error(['no metadata file found in ' dataDir]);
            else
                metafilename = fullfile(dataDir, listing(1).name);
                % exclude notes.txt file
                if strfind(metafilename,'notes')
                    metafilename = fullfile(dataDir, listing(2).name);
                end
                disp(['metadata file: ' metafilename]);
            end 

            % get the image filename format
            if ~exist(dataDir, 'dir')
                error(['data dir does not exist ' dataDir]);
            end
            
            [~,baremetafilename,~] = fileparts(metafilename);
            listing = dir(fullfile(dataDir,[baremetafilename '*.tif']));
            i = 1;

            if ~isempty(listing)
                filename = listing(i).name;
                while ( strcmp(filename(1),'.') || ~isempty(strfind(filename,'Preview')) )...
                        && i < numel(listing)
                    i = i + 1;
                    filename = listing(i).name;
                end
                [s,matches] = strsplit(filename,{'_.[0-9]{4}'},...
                                       'DelimiterType','RegularExpression',...
                                       'CollapseDelimiters',false);
                for i = 1:numel(matches)
                    matches{i} = [matches{i}(1:2) '%.4d'];
                end
                if numel(s) > 1
                    this.filename = [s{1} matches{:} s{end}];
                else
                    this.filename = filename;
                end
            else
                warning('Could not find image files');
            end

            % open the meta data file           
            fid = fopen(metafilename);

            if fid == -1
                error('file not found');
            end

            % move to the first line of the heads
            for i = 1:3
                tline = fgets(fid);
            end

            % scan through the header
            while ischar(tline) 
                k = strfind(tline,' : ');
                if isempty(k)
                    break;
                else
                    val = strtrim(tline(k + 3:end));

                    if ~isnan(str2double(val))
                        rawMeta.(tline(1:k-1)) = str2double(val);
                    else
                        rawMeta.(tline(1:k-1)) = val;
                    end
                end
                tline = fgets(fid);
            end

            % first make a field tPerFile to keep file information together 
            if isfield(rawMeta,'Time')
                this.tPerFile = NaN;
            end

            % resolution
            k = strfind(rawMeta.x,'*');
            this.xres = str2double(rawMeta.x(k+2:k+9));
            this.yres = str2double(rawMeta.y(k+2:k+9));

            % read out number of slices
            nZslices = 1;
            if isfield(rawMeta,'Z')
                nZslices = rawMeta.Z;
            end
            this.nZslices = nZslices;

            % read out number of channels
            nChannels = 1;
            if isfield(rawMeta,'Wavelength')
                nChannels = rawMeta.Wavelength;
            end
            this.nChannels = nChannels;

            % read out number of time points
            nTimePts = 1;
            if isfield(rawMeta,'Time')
                nTimePts = rawMeta.Time;
            end
            this.nTime = nTimePts;

            % read out number of positions
            nPositions = 1;
            if isfield(rawMeta,'Montage')
                nPositions = rawMeta.Montage;
            elseif isfield(rawMeta,'XY')
                nPositions = rawMeta.XY;
            end
            this.nPositions = nPositions;

            % read time interval
            if nTimePts > 1 
                while ischar(tline) 
                    tline = fgets(fid);
                    k = strfind(tline, 'Repeat T');
                    if ~isempty(k)
                        break;
                    end
                end
                s = strsplit(tline,{'(',')'});
                this.timeInterval = s{2};
            end

            % montage grid size
            if nPositions > 1 && isfield(rawMeta,'Montage')
                while ischar(tline) 
                    tline = fgets(fid);
                    k = strfind(tline, 'Montage Positions');
                    if ~isempty(k)
                        break;
                    end
                end
                s = strsplit(tline,' ');
                si = find(strcmp(s,'by'));
                this.montageGridSize = [str2double(s{si-1}) str2double(s{si+1})];
            end

            % channel names
            while ischar(tline) 
                tline = fgets(fid);
                k = strfind(tline, 'Repeat - Channel');
                if ~isempty(k)
                    break;
                end
            end
            s = strsplit(tline,{'(',')'});
            s = strsplit(s{2},',');
            this.channelNames = s;
            
            % read positions        
            if nPositions > 1
                
                if isfield(rawMeta,'Montage')        
                    
                    % overlap for montage
                    while ischar(tline) 
                        tline = fgets(fid);
                        k = strfind(tline, 'Overlap');
                        if ~isempty(k)
                            break;
                        end
                    end
                    overlap = str2double(strtrim(tline(k+8:end)));
                    this.montageOverlap = overlap;
                end

                % move on to XYZScan to read positions
                while ischar(tline) 
                    if ~isempty(strfind(tline,'[Region Info (Fields)]'))
                        break;
                    end
                    tline = fgets(fid);
                end
                for i = 1:3
                    tline = fgets(fid);
                end

                % read positions
                XYZ = zeros([nPositions 3]);
                for i = 1:nPositions
                    s = strsplit(tline,'\t');
                    XYZ(i,1) = str2double(s{end-1});
                    XYZ(i,2) = str2double(s{end});
                    tline = fgets(fid);
                    s = strsplit(tline,'\t');
                    XYZ(i,3) = str2double(s{1});
                    tline = fgets(fid);
                end
                this.XYZ = XYZ;
            end

            fclose(fid);

            % number of time points per file
            %---------------------------------
            if isfield(rawMeta,'Time') && exist('filename','var')

                pathstr = dataDir;%fileparts(filename);

                % determine the number of files the time series is split up into
                listing = dir(fullfile(pathstr,'*tif'));
                info = [];
                fileTmax = 0;
                for i = 1:numel(listing)
                    k = strfind(listing(i).name,'_t');
                    if ~isempty(k) && ~strcmp(listing(i).name(1),'.')
                        fileTmax = max(fileTmax, str2double(listing(i).name(k+2:k+5)));
                        if isempty(info)
                            info = imfinfo(fullfile(pathstr,listing(i).name));
                        end
                    end
                end
                
                % if files don't have a _t label nothing was found
                if isempty(info)
                    this.tPerFile = this.nTime;
                else
                    % time points per file
                    this.tPerFile = numel(info)/this.nZslices;
                    % for the case where wavelenghts were not exported
                    % separately
                    if isempty(strfind(this.filename,'_w'))
                        this.tPerFile = this.tPerFile/this.nChannels;
                    end
                end
            else
                this.tPerFile = 1;
            end
            
            this.raw = rawMeta;
        end
    end
end