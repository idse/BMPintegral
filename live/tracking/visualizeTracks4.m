function visualizeTracks4(positions, opts)
%Track input is a tree object holding compound tracks with lineages
%Extract cell data from position object

points = transpose({positions.cellData.XY});
if ~opts.use_tif
    areas = transpose({positions.cellData.area});
end

nucChannel = positions.nucChannel + 1;
trees = positions.tree;

% parse input options
if ~isfield(opts,'write_to_avi')
    opts.write_to_avi = true;
end
if ~isfield(opts,'write_to_tif')
    opts.write_to_tif = false;
end
if ~isfield(opts,'visible')
    opts.visible = 'off';
end
if ~isfield(opts,'default_tif')
    opts.default_tif = false;
end

%read the tif (MIP) file
if opts.use_tif
    if opts.default_tif
        path = opts.dataDir;
        name = positions.filename;
        name = strcat(name,'000',num2str(nucChannel-1));
        fname = fullfile(path, strcat(name,'.tif'));
    else
        [selectfile, selectpath] = uigetfile('*.tif','Select tif file');
        fname = fullfile(strcat(selectpath, selectfile));
        [path, name, ~] = fileparts(fname);
    end
    t = imfinfo(fname);
    %is there a way to get imsize without reading the tif?
    imsize = [t(1).Height, t(1).Width];
else
    max_pos = ceil(max(cellfun(@(x) max(x,[],'all'), points)));
    imsize = [max_pos, max_pos];
    path = uigetdir(pwd, 'Select MIP folder');
end

%load the segmentation
if opts.use_seg && opts.use_tif    
    if opts.default_seg
        segname = fullfile(path,[name,'_Simple Segmentation.h5']);
    else
        [segfile, segpath] = uigetfile('*.h5','Select segmentation file');
        segname = fullfile(segpath, segfile);
    end
    seg = h5read(segname,'/exported_data');
    seg = squeeze(seg);
    if ~isempty(nucChannel)
        seg = seg == nucChannel;
    else
        seg = seg == 1;
    end
    if ndims(seg) == 3
        seg = permute(seg,[2,1,3]);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %clean up the segmentation/use clean segmentation?
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%Find the number of compound tracks as the number of children of the root
%of the tree
track_nodes = getchildren(trees,1);
ntracks = length(track_nodes);
ntimepoints = size(points,1);

%Define an array of colors for visualization
%Distuinguishable colors function makes at most 9000 different colors; if
%there are > 9000 tracks, start repeating colors
if ntracks <= 9000
    colors = distinguishable_colors(ntracks,{'g','k'});
else
    colors = distinguishable_colors(9000);
    numreps = ceil(ntracks/9000) - 1;
    remainder = mod(ntracks,9000); % = ntracks - numreps*9000
    colors = vertcat(repmat(colors,numreps), colors(1:remainder,:));
end

%Pick a filename if writing to a tif file (or use the default)
if opts.write_to_tif || opts.write_to_avi
    writepath = fullfile(path,'Tracks');
    if ~exist(writepath,'dir')
        mkdir(writepath);
    end
    if opts.default_name && opts.use_tif
        filename = fullfile(writepath,[name, '_tracked_', date, '.tif']);
    else
        if isfield(opts,'name')
            filename = opts.name;
        else
            filename = input("Enter filename: ",'s');
        end
        %if a .tif extension is not present, add it (also remove any other
        %extension)
        if ~all(ismember('.tif', filename))
            dot_idx = find(filename == '.',1);
            if ~isempty(dot_idx)
                filename = strcat(filename(1:dot_idx-1),'.tif');
            else
                filename = strcat(filename,'.tif');
            end
        end
        filename = fullfile(writepath, filename);
    end
    aviname = strcat(filename(1:end-4),'.avi');
    
    if ~isfield(opts,'bigtiff')
        if ntimepoints > 600
            opts.bigtiff = true;
        else
            opts.bigtiff = false;
        end
    end
end
%if writing as an avi video, create a VideoWriter object
if opts.write_to_avi
    v = VideoWriter(aviname);
    open(v)
end

%set tags for the saved tif
if opts.write_to_tif
    tagstruct.Photometric = Tiff.Photometric.RGB;
    tagstruct.Compression = Tiff.Compression.None;
    tagstruct.BitsPerSample = 8;
    tagstruct.SamplesPerPixel = 3;
    tagstruct.ImageLength = imsize(1);
    tagstruct.ImageWidth = imsize(2);
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
end

%Make a cell array for visualizing the tracks that will hold xy positions
%versus time
tracks = cell(ntracks, 1);
%Parse the input tree to put it in a better format for visualization
for i = 1:ntracks
    node = track_nodes(i);
    subtree = buildsubtree(trees, node);
    tracks{i} = NaN(ntimepoints, 2*nnodes(subtree));
    for j = 1:nnodes(subtree)
        tracks{i}(:, 2*j-1) = get(subtree, j);
        for k = 1:ntimepoints
            idx = tracks{i}(k, 2*j-1);
            if ~isnan(idx)
                xy = points{k}(idx,:);
                tracks{i}(k, [2*j-1,2*j]) = xy;
            end
        end
    end
end

%% Do the visualization
f1 = figure('Units', 'normalized', 'Outerposition', [0, 0.04, 1, 0.96],...
    'Visible', opts.visible);
for j = 1:ntimepoints
    if mod(j,ceil(ntimepoints/20)) == 0
        disp(strcat(num2str(round(j/ntimepoints*100)),"%"))
    end
    clf %clear the current figure contents
    if opts.use_tif
        frame = double(imadjust(imread(fname,j))); %read the image at frame j, if using a tif
        frame = frame/max(frame,[],'all');
    else
        frame = ones(imsize(1), imsize(2)); %otherwise use a white background
    end
    ax = axes('Parent',f1);
    
    if opts.use_seg && opts.use_tif
        %if using the cell segmentation, overlay it in green on the image
        imz = cat(3, frame, frame + 0.3*seg(:,:,j), frame);
    else
        %otherwise, make the image rgb, but with no actual colors (allows a
        %black and white image with color tracks on top)
        imz = cat(3, frame, frame, frame);
    end
    
    imagesc(ax, imz);
    ax.PlotBoxAspectRatio = [1 1 1];
    if ~opts.use_tif
        %if not using a background image, visualize cells as circles
        radii = sqrt(areas{j}/pi);
        viscircles(points{j}, radii);
    elseif opts.use_tif && ~opts.use_seg
        %if not using a segmentation, mark the middle of cells with a small
        %circle/dot
%         viscircles(points{j}, 3*ones(size(points{j},1),1));
        hold(ax, 'on')
        scatter(points{j}(:,1),points{j}(:,2),10,'red','filled');
    end
    
    for i = 1:length(tracks)
        track = tracks{i};
        for k = 1:(size(track,2)/2)
            sidx = max(1, j - opts.tail_length);
            x = track(sidx:j,2*k-1);
            x = x(~isnan(x));
            y = track(sidx:j,2*k);
            y = y(~isnan(y));
            if ~isempty(x)
                line(ax,x,y,'Color',colors(i,:),'LineWidth',1)
            end
        end
    end
    
    if opts.write_to_tif || opts.write_to_avi
        frame = getframe(ax);
        frame = frame.cdata;
    end
    
    %write result to a tiff file
    if opts.write_to_tif
        if j == 1
            if opts.bigtiff
                disp('Writing as Big Tiff')
                mode = 'w8';
            else
                mode = 'w';
            end
        else
            mode = 'a';
        end
        t = Tiff(filename,mode);
        tagstruct.ImageLength = size(frame,1);
        tagstruct.ImageWidth = size(frame,2);
        t.setTag(tagstruct);
        t.write(frame);
        close(t)
    end
    
    if opts.write_to_avi
        writeVideo(v, frame);
        if j == ntimepoints
            %Close the VideoWriter after the last frame
            close(v);
        end
    end
    if strcmp(opts.visible,'on')
        title(strcat("Time = ", num2str(j)))
        pause(0.05)
    end
end

end