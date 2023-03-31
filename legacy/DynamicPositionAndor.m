classdef DynamicPositionAndor < Position
    % Data class to store cell data in a field of view

    % ---------------------
    % Idse Heemskerk, 2016
    % ---------------------

    properties

        tPerFile
    end
    
    methods

        function this = DynamicPositionAndor(meta, ID)
            % constructor 
            %
            % DynamicPositionAndor(meta, ID)
            %
            % meta:     MetadataAndor object
            % ID:       position index
            
            % matlab sucks
            if nargin == 0
                return
            end
            
            this.nChannels = meta.nChannels;
            this.cellData = struct();
            this.nucChannel = meta.nucChannel;
            
            this.nTime = meta.nTime;
            if isempty(meta.tPerFile) %~isfield(struct(meta),'tPerFile') || 
                this.tPerFile = this.nTime;
                warning('tPerFile not in meta, assuming tPerFile=nTime');
            else
                this.tPerFile = meta.tPerFile;
            end
            this.ID = ID;
            
            % this is a clunky way to only put the actual position in the
            % filename format but leave the rest of the %.4d pieces
            filenameFormat = meta.filename;
            barefname = sprintf(filenameFormat,this.ID-1);
            barefname = barefname(1:end-2);
            filenameFormat(1:numel(barefname)) = barefname;
            this.filenameFormat = filenameFormat;
        end

        % saving and loading
        %---------------------------------

        function img = loadImage(this, dataDir, channel, time, useMIP)
            % load image of colony
            %
            % img = loadImage(dataDir, channel, time)
            %
            % dataDir:  main data directory 
            % channel:  desired channel to be loaded
            %
            % img:      loaded image
            %
            % for now assume Andor format input, can be expanded later

            % fti : time index of file, e.g. if tPerFile = 2
            % Andor times start at 0, our time index and subti at 1
            % subti : time index within file
            %
            % so time=1: fti = 0, subti = 1
            % time=30: fti = 16, subti = 1
            
            if ~exist('useMIP','var')
                useMIP = false;
            end
            
            if useMIP
                error('useMIP not implemented for DynamicPositionAndor');
            end
            
            warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffErrorAsWarning');

            nFiles = ceil(this.nTime/this.tPerFile);
            fti = ceil(time/this.tPerFile) - 1;
            subti = rem(time-1,this.tPerFile) + 1;
            
            % if channels were separated
            if strfind(this.filenameFormat,'_w')

                % there is no _t part if all timepoints fit in a single file
                if this.nTime > this.tPerFile
                    fname = fullfile(dataDir, sprintf(this.filenameFormat,fti,channel));
                else
                    fname = fullfile(dataDir, sprintf(this.filenameFormat,channel));
                end

                info = imfinfo(fname);
                w = info.Width;
                h = info.Height;
                if fti + 1 < nFiles || this.nTime==this.tPerFile
                    nZslices = numel(info)/this.tPerFile;
                else
                    nZslices = numel(info)/rem(this.nTime,this.tPerFile);
                end

                ioffset = (subti-1)*nZslices;
                img = zeros([h w nZslices],'uint16');
                for i = 1:nZslices
                    img(:,:,i) = imread(fname, ioffset + i);
                end
                
            % if channels were not separated
            else
                
                % there is no _t part if all timepoints fit in a single file
                if this.nTime > this.tPerFile
                    fname = fullfile(dataDir, sprintf(this.filenameFormat,fti));
                else
                    fname = fullfile(dataDir, sprintf(this.filenameFormat));
                end
                
                info = imfinfo(fname);
                w = info.Width;
                h = info.Height;
                if fti + 1 < nFiles || this.nTime==this.tPerFile
                    nZslices = numel(info)/(this.nChannels*this.tPerFile);
                else
                    nZslices = numel(info)/(this.nChannels*rem(this.nTime,this.tPerFile));
                end

                % order in unseparated Andor files is czt
                ioffset = (subti-1)*nZslices*this.nChannels;
                img = zeros([h w nZslices],'uint16');
                for i = 1:nZslices
                    img(:,:,i) = imread(fname, ioffset + this.nChannels*(i-1) + channel);
                end
            end
            
            %disp(['loaded image ' fname]);
        end
    end
end