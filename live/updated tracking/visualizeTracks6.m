function visualizeTracks6(positions, opts)
%Track input is a tree object holding compound tracks with lineages
%Extract cell data from position object

points = transpose({positions.cellData.XY});

Div = transpose({positions.cellData.divs});

nucChannel = positions.nucChannel + 1;
trees = positions.tree;

% parse input options
if ~isfield(opts,'write_to_avi')
    opts.write_to_avi = true;
end
if ~isfield(opts,'visible')
    opts.visible = 'off';
end

%set name (for possible avi saving) and figure out image size
writepath = fullfile(opts.dataDir,'Tracks');
if ~exist(writepath,'dir')
        mkdir(writepath);
end
aviname = fullfile(writepath,[positions.filename,...
    '000',num2str(nucChannel-1),'_tracked_', date, '.avi']);
%if writing as an avi video, create a VideoWriter object
if opts.write_to_avi
    disp(strcat("Writing results as video to ", aviname))
    v = VideoWriter(aviname);
    open(v)
end

%set imsize (is this needed?)
testim = loadImage(positions,opts.dataDir,nucChannel-1,1);
imsize = [size(testim,1), size(testim,2)];

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

if isfield(opts,'xl')
    xl = opts.xl;
    yl = opts.yl;
else
    xl = [0.5, imsize(1) + 0.5];
    yl = [0.5, imsize(2) + 0.5];
end

%% Do the visualization
f1 = figure('Units', 'normalized', 'Outerposition', [0, 0.04, 1, 0.96],...
    'Visible', opts.visible);
for j = 1:ntimepoints
    divs = Div{j};
    if mod(j,ceil(ntimepoints/20)) == 0
        fprintf('.')
    end
    clf %clear the current figure contents
    frame = loadImage(positions,opts.dataDir,nucChannel-1,j);
    frame = im2double(imadjust(max(frame,[],3)));
    imz = cat(3, frame, frame, frame);
    
    ax = axes('Parent',f1);
    imagesc(ax, imz);
    ax.PlotBoxAspectRatio = [1 1 1];
    
    %mark the middle of cells with a small circle/dot
    hold(ax, 'on')
    scatter(points{j}(divs == 0,1),points{j}(divs == 0,2),10,'g','filled')
    scatter(points{j}(divs == 1,1),points{j}(divs == 1,2),10,'r','filled')
    %show the tails for tracks
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
    xlim(xl)
    ylim(yl)
    
    if opts.write_to_avi
        frame = getframe(ax);
        frame = frame.cdata;
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
fprintf('\n')

end