function trees = create_tracks_v2(this, opts)

%%%% Inputs %%%%
% this                      Position object with cell data
% opts                      Structure array with following fields:
%     max_linking_distance      maximum distance for frame-frame cell 
%                               movement
%     max_gap_distance          maximum cell movement distance when 
%                               doing gap closing (distance moved over 
%                               multiple frames)
%     max_gap_time           maximum # of frames for gap-closing, 
%                               merging, and splitting
%     imsize                    size of input image (get this from 
%                               meta instead?)

max_linking_distance = opts.max_linking_distance;
max_gap_distance = opts.max_gap_distance;
max_gap_time = opts.max_gap_time;
imsize = opts.imsize;
if isfield(opts,'tmax') && ~isempty(opts.tmax)
    tmax = opts.tmax;
else
    tmax = this.nTime;
end
if ~isfield(opts,'validate_links')
    opts.validate_links = false;
end
if ~isfield(opts,'dejitter')
    opts.dejitter = false;
end

%%%% Extract data from this %%%%
if ~opts.dejitter
    points = transpose({this.cellData.XY});
else
    disp('Dejittering')
    [points, ~] = dejitter(this, opts.dataDir);
    %save existing track position data to return to cellData at the end of
    %tracking
    oldXY = transpose({this.cellData.XY});
    for ti = 1:tmax
        this.cellData(ti).XY = points{ti};
    end
end

if isfield(this.cellData,'nucArea')
    areas = transpose({this.cellData.nucArea});
else
    %older version of Positions had area field instead of nucArea
    areas = transpose({this.cellData.area});
end
%identify dividing cells
for ti = 1:tmax
    %use object classification from Ilastik to identify dividing cells
    if isfield(this.cellData,'labels')
        labels = this.cellData(ti).labels;
        this.cellData(ti).divs = labels == 3;
    %if no Ilastik segmentation, but cell geometry statistics have been
    %computed, use these to estimate which cells are dividing
    elseif isfield(this.cellData,'nucMajorAxis')
        majors = this.cellData(ti).nucMajorAxis;
        minors = this.cellData(ti).nucMinorAxis;
        area = this.cellData(ti).nucArea;
        axisRatios = majors./minors;
        %dividing nuclei appear elongated, small, and bright (but nuclear
        %intensity is highly variable between cells, so don't use here)
        divs = axisRatios >= 2.5;
        divs = divs & (area < mean(area) - 2*std(area));
        this.cellData(ti).divs = divs;
    %otherwise treat all nuclei as not dividing
    else
        this.cellData(ti).divs = zeros(this.ncells(ti),1);
    end
end
Divs = transpose({this.cellData.divs});

%get the angle of the major axis of each cell wrt the right horizontal
%axis (i.e. positive x-axis) in radians
angles = cell(size(areas));
if isfield(this.cellData,'nucOrientation')
    for ti = 1:tmax
        angles{ti} = this.cellData(ti).nucOrientation*pi/180 + pi/2;
    end
else
    for ti =1:tmax
        angles{ti} = zeros(this.ncells(ti),1);
    end
end

%get nuclear intensity for each cell in each frame
intensities = transpose({this.cellData.nucLevel});
for i=1:length(intensities)
    intensities{i} = intensities{i}(:,this.nucChannel + 1);
end
ntimepoints = tmax;

%initialize simple tracks
% simple_tracks = cell(size(points{1},1),1);
% for i=1:size(points{1},1)
%     simple_tracks{i} = [i, 1];
% end
%this is faster but much less readable -> initialize each track with one
%point index according to # of cells in first frame, and time index of 1
simple_tracks =...
    num2cell([(1:size(points{1},1))',ones(size(points{1},1),1)],2);
%goods is used to keep track of which tracks are still active; for each
%frame, if a track i ended (linked to -1), then goods(i) is set to zero,
%and if a new track is started, a new array is added to simple_tracks and a
%new value of 1 is appended to goods
goods = ones(size(simple_tracks));

%%%% Frame-frame linking %%%%
disp('Doing frame-frame linking step')
tic
%Variables to hold data during linking
target_indices = cell(size(points));
source_indices = cell(size(points));
%Frame-to-frame linking step
f = figure;
for ti=1:ntimepoints - 1
    [ target_idxs, source_idxs ] = ...
        framelinker(this, ti, imsize, max_linking_distance);
    target_indices{ti} = target_idxs;
    source_indices{ti} = source_idxs;
    if opts.validate_links
        if ti == 1
            im1 = loadImage(this, opts.dataDir, this.nucChannel, ti);
            im1 = imadjust(max(im1,[],3));
        else
            im1 = im2;
        end
        im2 = loadImage(this, opts.dataDir, this.nucChannel, ti+1);
        im2 = imadjust(max(im2,[],3));
        cla
        imshowpair(im1,im2);
        if ti == 1
            bounds = f.OuterPosition;
            bounds(1) = 0;
            f.OuterPosition = bounds;
        end
        hold on
        x = points{ti}(:,1);
        y = points{ti}(:,2);
        scatter(x, y, 36, 'g', 'filled')
        scatter(points{ti+1}(:,1),points{ti+1}(:,2),36,'b','filled')
        links = find(target_indices{ti} > 0);
        nonlinks = find(target_indices{ti} == -1);
        u1 = points{ti+1}(target_indices{ti}(links),1) - x(links);
        v1 = points{ti+1}(target_indices{ti}(links),2) - y(links);
        u2 = zeros(length(nonlinks),1);
        v2 = zeros(length(nonlinks),1);
        quiver([x(links);x(nonlinks)],[y(links);y(nonlinks)],...
            [u1;u2],[v1;v2],0,'r','LineWidth',2)
        scatter(x(Divs{ti} == 1),y(Divs{ti} == 1),36,'r','filled')
        legend('source','target','links','dividing')
        title(num2str(ti))
        [lindx, ~] = listdlg('ListString',{'Continue','Break'});
        if lindx == 2
            opts.validate_links = false;
        end
    end
    %String tracks together as frames are being linked (doing this
    %keeps us from being able to use a parfor loop here though)
    %In each frame, only look at tracks for which goods(j) = 1; this gives
    %the indices of all such tracks in simple_tracks
    source_tracks = find(goods == 1);
    for j=1:length(source_tracks)
        jj = source_tracks(j); %index of the source track
        sp = simple_tracks{jj}(end,1); %cell index at last time in source track
        tp = target_indices{ti}(sp); %target cell for this cell
        %If the point in track j is linked to nothing in the next frame,
        %set goods(jj) = 0; this means the track ends in this frame
        if tp == -1
            goods(jj) = 0;
        else
            %Otherwise add the index of the point in the next frame to the
            %row of simple_tracks(jj) alongside its time index
            simple_tracks{jj} = [simple_tracks{jj};[tp, ti+1]];
        end
    end
    %Find points at time i which are linked to nothing in the previous
    %frame; these correspond to new tracks
    new_idxs = find(source_indices{ti} == -1);
    for j=1:length(new_idxs)
        jj = new_idxs(j);
        %Add an entry to simple tracks starting with the index and time of
        %the first point in the track
        simple_tracks{size(simple_tracks,1)+1} = [jj, ti+1];
        %Append a value of 1 to goods to indicate that this track is active
        goods = [goods; 1];
    end

    if mod(ti,round((ntimepoints-1)/20)) == 0
        fprintf('.')
    end
end
fprintf('\n')
toc
close(f)

ntracks = size(simple_tracks, 1);
disp(strcat("Number of simple tracks = ", num2str(ntracks)))

%suppress warning about temporary variables in parfor loop:
warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary')

%%%% Gap closing, splitting, merging, creating compound tracks %%%%
%this function should probably just be made to take the whole position
%object as an input and use that to extract relevant features within
%the function instead of passing 1000000 arguments
[gaps, merges, splits, ~, ~, ~] = newtracklinker(simple_tracks,...
     max_gap_time, max_gap_distance, points, intensities,...
     areas, Divs, angles);
 %turn the warning back on at the end
warning('on', 'MATLAB:mir_warning_maybe_uninitialized_temporary')
 %Output: nx4 arrays with linking info
%gaps(i,:) = [end_track, start_track, gap_time, cost]
%merges(i,:) = [end_track, mid_track, merge_time, cost]
%splits(i,:) = [end_track, mid_track, split_time, cost]


%%%% Close gaps in simple_tracks %%%%
disp("Closing gaps in simple tracks")
tCol = 2;
idxCol = 1;
%Start by grouping gap-closed tracks together - track_idxs should hold a
%cell array where indices of these tracks are grouped; i.e. if track 3 is
%closed to track 12 and the end of track 12 links to the beginning of track
%30 and track 30 terminates without another link, then one of the entries
%of track_idxs will be [3,12,30], and the track_ids for 12 and 30 will be
%changed to 0
%track_idxs{i,2} holds a nx2 array with costs and times for gap
%closes done between the simple_tracks in track_idxs{i,1}
track_idxs = cell(ntracks,1);
gap_conf = cell(ntracks,1);
track_ids = ones(ntracks,1);
for i=1:ntracks
    track_ids(i) = 1; %unnnecessary?
    track_idxs{i,1} = i;
end

for i = 1:ntracks
    idx = i;
    if ismember(idx,gaps(:,1)) && track_ids(idx)
        idxlist = idx;
        cts = []; %costs and times
        criterion = true;
        while criterion == true
            idxNext = gaps(gaps(:,1) == idx,2);
            cost_time = gaps(gaps(:,1) == idx,3:4);
            track_ids(idxNext) = 0;
            idxlist = [idxlist,idxNext];
            cts = [cts; cost_time];
            if ismember(idxNext,gaps(:,1))
                if idx == idxNext
                %prevent infinite while loop when there is an error in CM construction
                    error('Incorrect track assignments')
                end
                idx = idxNext;
            else
                criterion = false;
            end
        end
        track_idxs{i,1} = idxlist;
        gap_conf{i,1} = cts;
    end
end

%Revise the track_idxs to get rid of duplicates using track_ids and reset
%track_ids for merging and splitting
track_idxs2 = cell(sum(track_ids),1);
gap_conf2 = cell(sum(track_ids), 1);
idx = 1;
for i=1:size(track_idxs,1)
    if track_ids(i) == 1
        track_idxs2{idx,1} = track_idxs{i,1}; %track indices
        gap_conf2{idx,1} = gap_conf{i,1}; %costs for linking & times
        idx = idx + 1;
    end
end
track_idxs = track_idxs2;
gap_conf = gap_conf2;
track_ids = ones(size(track_idxs,1),1);

%Make gap_tracks using track_idxs; this makes groups of tracks that are
%gap-closed into single tracks with discontinuities
gap_tracks = cell(length(track_ids),1);
for i=1:size(gap_tracks,1)
    %time points with no cell have nan for the cell index
    gap_tracks{i} = nan(ntimepoints,1);
    for j=1:length(track_idxs{i,1})
        sidx = track_idxs{i,1}(j); %index of corresponding simple track
        times = simple_tracks{sidx}(:,tCol); %times for which the track is defined
        idxs = simple_tracks{sidx}(:,idxCol); %cell index in frame at these times
        gap_tracks{i}(times) = idxs;
    end
end
%gap_conf{i,1} still has the costs and times for gap_tracks{i}

%Go through merges and splits and replace track indices of simple_tracks
%with corresponding track indices for gap_tracks
splits2 = splits;
merges2 = merges;
for i = 1:size(splits2,1)
    for j = 1:2
        sidx = splits2(i,j);
        gidx = find(cellfun(@(x) ismember(sidx,x), track_idxs));
        splits2(i,j) = gidx;
    end
end

for i = 1:size(merges2,1)
    for j = 1:2
        sidx = merges2(i,j);
        gidx = find(cellfun(@(x) ismember(sidx,x), track_idxs));
        merges2(i,j) = gidx;
    end
end

%%%% Handle merging events %%%%
%Merging events are assumed to be artifacts of the detection, segmentation,
%or tracking process and should be handled before using splitting events to
%create compound tracks

disp("Resolving merging and splitting events")
%Find which tracks are merged to
merge_tracks = unique(merges2(:,2));

%For each track that is merged to, create a list of events that happen to
%it: merges and splits, the tracks that are merging and splitting to/from
%them, and the time point at which event occurs; label merges and splits
%with 0 and 1, respectively
%[merge_track, other track, time, merge/split (0/1)]
events = cell(length(merge_tracks),1);
for k = 1:length(merge_tracks)
    mevents = merges2(merges2(:,2) == merge_tracks(k),[2,1,3]);
    sevents = splits2(splits2(:,2) == merge_tracks(k),[2,1,3]);
    events{k} = [mevents, zeros(size(mevents,1),1); sevents, ones(size(sevents,1),1)];
    [~, I] = sort(events{k}(:,3));
    events{k} = events{k}(I,:);
end

%For each unique track that experiences at least one merging event, create
%a cost matrix to make assignments: each of the two branches in a merging
%event may be assigned either to a subsequent split branch or considered an
%error
%Set a cutoff such that if there are more than so many frames after a
%merge before another split, it is considereed an erroneous assignment

%For now, essentially just gap-close, assuming that if there is a split
%soon after a merge then the merged track corresponds to the split track
%and otherwise that the merged track terminates; this is not robust
merge_cutoff = 10;
[newtracks, splits4] = newmergelinker(events, gap_tracks,...
splits2, points, areas, intensities, merge_cutoff);

disp("Constructing final tracks")
trees = build_tree(newtracks, splits4);

nctracks = numel(getchildren(trees,1));
disp(strcat("Number of compound tracks = ", num2str(nctracks)))

%if cell location data was modified during dejittering, set it back
if opts.dejitter
    for ti = 1:tmax
        this.cellData(ti).XY = oldXY{ti};
    end
end




    
end