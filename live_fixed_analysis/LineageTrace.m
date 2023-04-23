classdef LineageTrace < handle
   
    properties
        live_position       %Position object from time-lapse data
        fixed_position      %Position object from fixed cell data
        
        liveMeta            %live image metadata
        fixedMeta           %fixed image metadata
        
        histories  %cell array with a struct for each position;
                            %each entry of the cell is a struct with fate
                            %marker values and signaling histories for
                            %cells in the fixed image
        
        %image size depends on stitching parameters, but all frames and
        %channels should have the same size in the live data and all
        %channels in the fixed data should have the same image size based
        %on the stitching approach (or [1024, 1024] if not stitched)
        liveSize            %Size (in pixels) of images in the time-lapse
        fixedSize           %Size (in pixels) of fixed images
        dataDir             %Base folder containing the data & script
                            %only used to help locate liveDir and fixedDir
        liveDir             %Folder containing raw data for time-lapse
        fixedDir            %Folder containing raw data for fixed images
        mapped_idxs         %Cell with list of indices of live-cell points
                            %mapped to cells in fixed data for each
                            %position; -1 indicates no corresponding cell
        shiftyx             %relative shift found between live and fixed
                            %images (cell with 1x2 entries)
    end
    
    methods
        
        %constructor
        function this = LineageTrace(varargin)
            %Possible constructor calls:
            % this = LineageTrace(live_dir, fixed_dir) %standard constructor
            % this = LineageTrace(positions_live, positions_fixed) %used if input class is 'Position'
            % this = LineageTrace(positions_live, positions_fixed, dataDir)
            % this = LineageTrace()
            
            %initialize with two position objects
            if nargin == 2
                if isa(varargin{1},'Position')
                    %LineageTrace(live_position, fixed_position)
                    this.live_position = varargin{1};
                    this.fixed_position = varargin{2};
                    this.dataDir = pwd;
                    this.selectLiveDir();
                    this.selectFixedDir();
                else
                    %LineageTrace(liveDir, fixedDir)
                    this.liveDir = varargin{1};
                    this.fixedDir = varargin{2};
                    S = load(fullfile(this.liveDir,'positions.mat'));
                    this.live_position = S.positions;
                    S = load(fullfile(this.fixedDir,'positions.mat'));
                    this.fixed_position = S.positions;
                end
            elseif nargin ==3
                %LineageTrace(live_position, fixed_position, dataDir)
                this.live_position = varargin{1};
                this.fixed_position = varargin{2};
                this.dataDir = varargin{3};
                this.selectLiveDir();
                this.selectFixedDir();
            elseif nargin == 4
                %LineageTrace(live_position, fixed_position, liveDir, fixedDir)
                this.live_position = varargin{1};
                this.fixed_position = varargin{2};
                this.liveDir = varargin{3};
                this.fixedDir = varargin{4};
            elseif nargin == 6
                %LineageTrace(live_position, fixed_position, liveMeta, fixedMeta, liveDir, fixedDir)
                this.live_position = varargin{1};
                this.fixed_position = varargin{2};
                this.liveMeta = varargin{3};
                this.fixedMeta = varargin{4};
                this.liveDir = varargin{5};
                this.fixedDir = varargin{6};
            %initialize with no arguments
            elseif nargin == 0
                %properties may be manually set
                return
            end
            
            if nargin > 0 && nargin < 6
                liveMeta = load(fullfile(this.liveDir,'meta.mat'));
                this.liveMeta = liveMeta.meta;
                fixedMeta = load(fullfile(this.fixedDir,'meta.mat'));
                this.fixedMeta = fixedMeta.meta;
            end
            
            this.histories = cell(size(this.live_position));
            this.mapped_idxs = [];
            this.shiftyx = cell(size(this.live_position));
        end
        
        %manually select live cell data folder
        function selectLiveDir(this)
            this.liveDir = uigetdir(this.dataDir, 'Select raw live data folder');
        end
        %manually select fixed data folder
        function selectFixedDir(this)
            this.fixedDir = uigetdir(this.dataDir, 'Select raw fixed data folder');
        end
        
        function XY = align_positions(this, pidx, alignmentMethod)
            %pi is the index of the position for which the images are being aligned
            %this is mostly called from within mapPoints instead of as a
            %standalone function
            %alignmentMethod specifies how images/datasets from subsequent
            %imaging rounds are aligned with the options:
                %automatic - use cross correlation to rigidly register
                %nuclear images
                %colonyCenter - for micropatterns, just align the colony
                %centroids
                %manual - if for some reason the automatic alignment
                %doesn't work well, there is the option to manually align
                %the images interactively
                %TODO: add the option to warp images onto one another to
                %account for colony deformation or cell movement
            %XY is a ncellsx2 array of adjusted cell locations in fixed image
            
            %if the fixed directory is not set, set it
            if isempty(this.fixedDir)
                selectFixedDir(this);
            end
            %if the live directory is not set, set it
            if isempty(this.liveDir)
                selectLiveDir(this);
            end
            
            p_live = this.live_position(pidx);
            p_fixed = this.fixed_position(pidx);
            %load images for image-based alignment
            try
                live = p_live.loadImage(this.liveDir, p_live.nucChannel, p_live.nTime);
            catch
                live = p_live.loadImage(fullfile(this.liveDir,'MIP'), p_live.nucChannel, p_live.nTime);
            end
            live = max(live,[],3); %take MIP in case of z-slices
            try
                fixed = p_fixed.loadImage(this.fixedDir, p_fixed.nucChannel, 1);
            catch
                fixed = p_fixed.loadImage(fullfile(this.fixedDir,'MIP'), p_fixed.nucChannel, 1);
            end
            fixed = max(fixed,[],3); %take MIP in case of z-slices
            %make images the same size
            this.liveSize = size(live);
            this.fixedSize = size(fixed);
            mn = min(size(live), size(fixed));
            m = mn(1); n= mn(2);
            live = live(1:m,1:n);
            fixed = fixed(1:m,1:n);
            
            if strcmp(alignmentMethod,'automatic')
                %automatically determine shift between fixed and live images
                [shiftx,shifty,~,~] = xcorr2fft(live, fixed);
            elseif strcmp(alignmentMethod,'colonyCenter')
                % alignment by colony center
                if ~isa(p_live,'Colony') || ~isa(p_fixed,'Colony')
                    error('colony centers can only be used for data of Colony type')
                end
                center1 = p_live.center;
                center2 = p_fixed.center;
                shift = round(center2 - center1);
                shifty = shift(1); shiftx = shift(2);
            elseif strcmp(alignmentMethod,'manual')
                %manually find the optimal image shift for each position
                liveIm = imadjust(live,stitchedlim(live));
                fixedIm = imadjust(fixed,stitchedlim(fixed));
                shiftx = 0;
                shifty = 0;
                close all
                f = figure('WindowState','maximized');
                p1 = imshow(cat(3,fixedIm,liveIm,fixedIm),'InitialMagnification','fit');
                cleanSubplot
                breakvar = false;
                while breakvar == false
                    xinrange = max(1,1+shiftx):min(n,n+shiftx);
                    xoutrange = max(1,1-shiftx):min(n,n-shiftx);
                    yinrange = max(1,1+shifty):min(m,m+shifty);
                    youtrange = max(1,1-shifty):min(m,m-shifty);

                    aligned = zeros(mn,'uint16');
%                     aligned(xoutrange, youtrange) = fixedIm(xinrange, yinrange);
                    aligned(youtrange, xoutrange) = fixedIm(yinrange, xinrange);
                    %update image for new frame
                    set(p1,'CData',cat(3,aligned,liveIm,aligned));

                    title(sprintf('shiftx = %d, shifty = %d', shiftx, shifty))

                    waitforbuttonpress;
                    key = f.CurrentCharacter;

                    if strcmp(key,'s')
                        shiftx = shiftx - 1;
                    elseif strcmp(key,'w')
                        shiftx = shiftx + 1;
                    elseif strcmp(key,'d')
                        shifty = shifty - 1;
                    elseif strcmp(key,'a')
                        shifty = shifty + 1;
                    elseif strcmp(key,'+')
                        xl = get(gca,'xlim');
                        yl = get(gca,'ylim');
                        set(gca,'xlim', mean(xl) + 0.25*(xl(2) - xl(1))*[-1 1])
                        set(gca,'ylim', mean(yl) + 0.25*(yl(2) - yl(1))*[-1 1])
                    elseif strcmp(key,'-')
                        xl = get(gca,'xlim');
                        yl = get(gca,'ylim');
                        set(gca,'xlim', mean(xl) + (xl(2) - xl(1))*[-1 1])
                        set(gca,'ylim', mean(yl) + (yl(2) - yl(1))*[-1 1])
                    elseif strcmp(key,'4')
                        shifty = shifty + 10;
                    elseif strcmp(key,'6')
                        shifty = shifty - 10;
                    elseif strcmp(key,'2')
                        shiftx = shiftx - 10;
                    elseif strcmp(key,'8')
                        shiftx = shiftx + 10;
                    elseif strcmp(key,'e')
                        breakvar = true;
                    end
                    drawnow limitrate
                end
                close(f)
                shiftx = -shiftx; shifty = -shifty;
            end
            this.shiftyx{pidx} = [shifty shiftx];
            %shift the XY positions of cells in the fixed frame accordingly
            fixed_XY = p_fixed.cellData.XY;
            XY = fixed_XY - [shifty, shiftx];
        end
        
        function mapPoints(this, max_distance, opts)
            %this.mapPoints(max_distance) 
            %this.mapPoints(max_distance, positionIdx) 
            %
            %maps the cells from the fixed
            %data to cells at the end of the live data for every position
            %max_distance is maximum expected cell movement distance in pixels
            %typical max_distance for 40x images = 30 (half for 20x)
            %if useZ, use xyz positions to determine the cost for linking,
            %with coordinates converted from pixels/slices to um; also need
            %to specify max_distance in um if using Z info
            
            if ~exist('opts','var')
                opts = struct;
            end
            
            if isfield(opts,'positionIdx')
                positionIdx = opts.positionIdx;
            else
                positionIdx = 1:numel(this.live_position);
            end
            
            if ~isfield(opts,'alignmentMethod')
                opts.alignmentMethod = 'automatic';
            end
            
            if isfield(opts,'useZ')
                useZ = opts.useZ;
            else
                useZ = false;
            end
            
            if useZ
                disp('Using Z coordinate and assuming max_distance is specified in microns')
            end
            
            this.mapped_idxs = cell(size(this.live_position));
            for pidx = positionIdx
                if useZ
%                     source_info = struct;
%                     target_info = struct;
%                     %centroid XY position
%                     source_info.XY = this.live_position(pidx).cellData(end).XY;
%                     target_info.XY =...
%                         align_positions(this, pidx, opts.alignmentMethod);
%                     %nuclear Z slice
%                     source_info.Z = this.live_position(pidx).cellData(end).nucZ;
%                     target_info.Z = this.fixed_position(pidx).cellData(1).nucZ;
%                     %nuclear area
%                     source_info.nucArea = this.live_position(pidx).cellData(end).nucArea;
%                     target_info.nucArea = this.fixed_position(pidx).cellData(1).nucArea;
%                     %metadata
%                     source_info.meta = this.liveMeta;
%                     target_info.meta = this.fixedMeta;
%                     
%                     imsize = min(this.liveSize, this.fixedSize);
%                     
%                     [ ~, source_indices, ~] =...
%                         LiveFixedLinkerWithZ(source_info, target_info, imsize, max_distance);
                    
                    %nucleus centroid xy positions
                    xy1 = this.live_position(pidx).cellData(end).XY;
                    xy2 = align_positions(this, pidx, opts.alignmentMethod);
                    %colony-centered and in microns
                    if isa(this.live_position,'Colony')
                        cm = this.live_position(pidx).center(end,:);
                    else
                        cm = [0 0];
                    end
                    
                    xy1 = this.liveMeta.xres*(xy1 - cm);
                    xy2 = this.fixedMeta.xres*(xy2 - cm);
                    %image size
                    imsize = this.liveMeta.xres*min(this.liveSize, this.fixedSize);
                    
                    Z1 = this.live_position(pidx).cellData(end).nucZ;
                    Z2 = this.fixed_position(pidx).cellData(1).nucZ;
                    u1 = mean(Z1); s1 = std(Z1);
                    z1 = this.liveMeta.zres*(Z1 - u1);
                    u2 = mean(Z2); s2 = std(Z2);
                    z2 = this.liveMeta.zres*s1/s2*(Z2 - u2);

                    sinfo = struct(...
                        'XY',       xy1,...
                        'Z',        z1,...
                        'nucArea',  this.live_position(pidx).cellData(end).nucArea);

                    tinfo = struct(...
                        'XY',       xy2,...
                        'Z',        z2,...
                        'nucArea',  this.fixed_position(pidx).cellData(1).nucArea);
                    
                    [~, source_indices, ~] =...
                        LiveFixedLinkerWithZ_v2(sinfo, tinfo, imsize, max_distance);
                    
                    this.mapped_idxs{pidx} = source_indices;
                    fprintf('Number of links = %d\n', sum(source_indices > 0))
                else
                    %get information for frame-frame cell linking
                    source_info = ones(this.live_position(pidx).ncells(end), 4); %initialize source_info
                    target_info = ones(this.fixed_position(pidx).ncells(1), 4); %initialize target_info
                    source_info(:,1:2) = this.live_position(pidx).cellData(end).XY;
                    
                    %align images to get cell positions for target_info
                    disp('Aligning positions')
                    target_info(:,1:2) =...
                        align_positions(this, pidx, opts.alignmentMethod);
                    imsize = min(this.liveSize, this.fixedSize);
                    
                    %do the linking as a linear assignment problem
                    disp('Linking cells from live to fixed data')
                    [ ~, source_indices, ~] = livefixedlinker(source_info,...
                        target_info, imsize, max_distance);
                    this.mapped_idxs{pidx} = source_indices;
                    fprintf('Number of links = %d\n', sum(source_indices > 0))
                end
            end
        end
        
        function checkAlignment(this, pidx)
            %after mapping points, overlay the (shifted) images and show
            %the link between fixed and live points
            p_live = this.live_position(pidx);
            p_fixed = this.fixed_position(pidx);
            
            try
                liveIm = p_live.loadImage(this.liveDir, p_live.nucChannel, p_live.nTime);
            catch
                liveIm = p_live.loadImage(fullfile(this.liveDir,'MIP'), p_live.nucChannel, p_live.nTime);
            end
            liveIm = max(liveIm,[],3);
%             liveIm = imadjust(liveIm,stitchedlim(liveIm));
            try
                fixedIm = p_fixed.loadImage(this.fixedDir, p_fixed.nucChannel, 1);
            catch
                fixedIm = p_fixed.loadImage(fullfile(this.fixedDir,'MIP'), p_fixed.nucChannel, 1);
            end
            fixedIm = max(fixedIm,[],3);
%             fixedIm = imadjust(fixedIm,stitchedlim(fixedIm));
            
            mn = min(size(liveIm), size(fixedIm));
            m = mn(1); n = mn(2);
            liveIm = liveIm(1:mn(1), 1:mn(2));
            fixedIm = fixedIm(1:mn(1), 1:mn(2));
            
            shifty = this.shiftyx{pidx}(1);
            shiftx = this.shiftyx{pidx}(2);
            
            xinrange = max(1,1+shiftx):min(m,m+shiftx);
            xoutrange = max(1,1-shiftx):min(m,m-shiftx);
            yinrange = max(1,1+shifty):min(n,n+shifty);
            youtrange = max(1,1-shifty):min(n,n-shifty);
            
            aligned = zeros(mn,'uint16');
            aligned(xoutrange, youtrange) = fixedIm(xinrange, yinrange);
            
            XYfixed = p_fixed.cellData(1).XY - [shifty, shiftx];
%             XYlive = p_live.cellData(p_live.nTime).XY;
            XYlive = p_live.cellData(end).XY;
            
            mappedSourceIdx = find(this.mapped_idxs{pidx} > 0);
            mappedTargetIdx = this.mapped_idxs{pidx}(this.mapped_idxs{pidx} > 0);
            dXY = XYlive(mappedTargetIdx,:) - XYfixed(mappedSourceIdx,:);
            
            imshowpair(imadjust(liveIm,stitchedlim(liveIm)),...
                imadjust(aligned,stitchedlim(aligned)))
            hold on
            scatter(XYfixed(:,1),XYfixed(:,2),10,'b','filled');
            scatter(XYlive(:,1),XYlive(:,2),10,'r','filled');
            quiver(XYfixed(mappedSourceIdx,1),XYfixed(mappedSourceIdx,2),dXY(:,1),dXY(:,2),0,'LineWidth',2);
            
        end
        
        function makeHistories(this, pidx, max_distance)
            %Align the nuclear image in the fixed data to the last frame of
            %the live data, then perform a nearest-neighbor linking to
            %establish correspondence between each cell (or most of them)
            %in the fixed data and a cell in the live data, and its
            %corresponding time trace
            
            if ~exist('max_distance','var')
                %default max distance should be made a function of the
                %image magnification
                max_distance = 30;
            end
            
            %specify fields from this.live_position.timeTraces that should
            %be saved in this.histories
            fields = {'XY','NCratio','nucLevel','cytLevel','nucArea','nucZ',...
                'nucMajorAxis','nucMinorAxis'};
            
            if isempty(this.mapped_idxs{pidx})
                this.mapPoints(max_distance, pidx)
            end
            
            fixedIdxs = find(this.mapped_idxs{pidx} > 0);
            liveIdxs = this.mapped_idxs{pidx}(fixedIdxs);
            
            G = this.live_position(pidx).G;
            %build tracks from digraph
            this.histories{pidx} = graphSignalingHistories(this.live_position(pidx), G, fields, liveIdxs);
            %assign fate markers
            fm = this.fixed_position(pidx).cellData.nucLevel(fixedIdxs,:);
            for ii = 1:length(this.histories{pidx})
                this.histories{pidx}(ii).fateMarkers = fm(ii,:);
            end
            
        end
        
        function M = histories2mats(this,positionIdx,opts)
            
            npos = length(this.live_position);
            if ~exist('positionIdx','var')
                positionIdx = 1:npos;
            end
            
            if ~exist('opts','var')
                opts = struct;
            end
            
            ntime = this.live_position(1).nTime;
            fields = {'XY','NCratio','nucLevel','cytLevel','nucArea','nucZ',...
                'nucMajorAxis','nucMinorAxis'};
            nfields = length(fields);
            
            M = struct;
            for fi = 1:nfields
                M.(fields{fi}) = cell(1,npos);
            end
            
            M.fateMarkers = cell(1,npos);
            M.positionIdx = cell(1,npos);
            
            for pidx = positionIdx
                nucChannel = this.live_position(pidx).nucChannel;
                opts.nucChannel = nucChannel;
                start_times = cellfun(@(x) x(1), {this.histories{pidx}.Time});
                hists = this.histories{pidx}(start_times==1);
                nhists = length(hists);
                
                for fi = 1:nfields
                    nc = size(hists(1).(fields{fi}),2);
                    X = NaN(ntime,nhists,nc);
                    for ci = 1:nc
                        for ii = 1:nhists
                            %cleanHistory_v2(histories,idx,opts)
                            opts.channel = ci-1;
                            opts.field = fields{fi};
                            X(:,ii,ci) =...
                                cleanHistory_v2(hists,ii,opts);
%                             X(:,ii,ci) =...
%                                 cleanHistory(hists,ii,ci-1,fields{fi},nucChannel);
                        end
                    end
                    M.(fields{fi}){pidx} = X;
                end
                
                nchannels = length(hists(1).fateMarkers);
                fm = NaN(nchannels,nhists);
                for ii = 1:nhists
                    fm(:,ii) = hists(ii).fateMarkers';
                end
                M.fateMarkers{pidx} = fm;
                M.positionIdx{pidx} = pidx*ones(1,nhists);
            end
            
            for fi = 1:length(fields)
                M.(fields{fi}) = cell2mat(M.(fields{fi}));
            end
            M.fateMarkers = cell2mat(M.fateMarkers);
            M.positionIdx = cell2mat(M.positionIdx);
            
        end
        
    end
    
    
    
    
    
    
end