function tracks = create_tracks(positions, max_linking_distance,...
    max_gap_distance, max_gap_closing, imsize)
%The function takes as its input a variable of type Position with field
%cellData containing ntimepoint-length arrays of data on x-y positions,
%nuclear intensities, and areas of cells at a position in the video
%P = positions(1);

%% Set parameters (make these inputs)
%Introduce some struct with optional input parameters for image size, max
%linking distance, max gap closing, and eventually options for the type of
%motion/approach for making cost functions, not just cutoffs

%% Extract data from positions
if ~isempty(positions.nucChannel)
    nucChannel = positions.nucChannel + 1;
else
    nucChannel = 1;
end

points = transpose({positions.cellData.XY});
areas = transpose({positions.cellData.area}); %this probably should not be unused
intensities = transpose({positions.cellData.nucLevel});
for i=1:length(intensities)
    intensities{i} = intensities{i}(:,nucChannel);
end
ntimepoints = size(points,1);

%get rid of empty time points (pad with NaNs)
emptyTimes = cellfun(@isempty, points);
if sum(emptyTimes == 1) > 0
    for i = find(emptyTimes)
        points{i} = [NaN NaN];
    end
end

%% Frame-frame linking
disp('Doing frame-frame linking step')
tic
%Variables to hold data during linking
target_indices = cell(size(points));
source_indices = cell(size(points));

%Frame-to-frame linking step
for i=1:ntimepoints - 1
    source = [points{i}, areas{i}, intensities{i}];
    target = [points{i+1}, areas{i+1}, intensities{i+1}];
    [ target_idxs, source_idxs ] = ...
        newhungarianlinker(source, target, imsize, max_linking_distance);
    target_indices{i} = target_idxs;
    source_indices{i} = source_idxs;
end
toc

%% String frame-to-frame links into tracks over the full time period
%Simple tracks is a cell object that holds tracks with no gap closing,
%linking or merging
simple_tracks = cell(size(points{1},1), 1);

%goods is used to keep track of which tracks are still active; for each
%frame, if a track i ended (linked to -1), then goods(i) is set to zero,
%and if a new track is started, a new array is added to simple_tracks and a
%new value of 1 is appended to goods
goods = ones(size(simple_tracks));
for i=1:size(points{1},1)
    simple_tracks{i} = [i, 1];
end

for i=1:ntimepoints - 1
    if mod(i,round((ntimepoints-1)/20)) == 0
        disp(strcat(num2str(round(i*100/(ntimepoints-1))),"%"))
    end
    %In each frame, only look at tracks for which goods(j) = 1; this gives
    %the indices of all such tracks in simple_tracks
    source_tracks = find(goods == 1);
    for j=1:length(source_tracks)
        jj = source_tracks(j);
        sp = simple_tracks{jj}(end,1);
        tp = target_indices{i}(sp);
        %If the point in track j is linked to nothing in the next frame,
        %set goods(jj) = 0; this means the track ends in this frame
        if tp == -1
            goods(jj) = 0;
        else
            %Otherwise add the index of the point in the next frame to the
            %row of simple_tracks(jj) alongside its time index
            simple_tracks{jj} = [simple_tracks{jj};[tp, i+1]];
        end
    end
    %Find points at time i which are linked to nothing in the previous
    %frame; these correspond to new tracks
    new_idxs = find(source_indices{i} == -1);
    for j=1:length(new_idxs)
        jj = new_idxs(j);
        %Add an entry to simple tracks starting with the index and time of
        %the first point in the track
        simple_tracks{size(simple_tracks,1)+1} = [jj, i+1];
        %Append a value of 1 to goods to indicate that this track is active
        goods = [goods; 1];
    end
end

ntracks = size(simple_tracks, 1);
disp(strcat("Number of simple tracks = ", num2str(ntracks)))

%% Gap closing, splitting, merging, creating compound tracks
[gaps, merges, splits, ~] = tracklinker5(simple_tracks,...
     max_gap_closing, max_gap_distance, points, intensities);
 
% [gaps, merges, splits, ~] = tracklinker3(simple_tracks,...
%      max_gap_closing, max_gap_distance, points, intensities);

%ntracks = size(simple_tracks, 1);

%% Close gaps in simple_tracks
tCol = 2;
idxCol = 1;

%Start by grouping gap-closed tracks together - track_idxs should hold a
%cell array where indices of these tracks are grouped; i.e. if track 3 is
%closed to track 12 and the end of track 12 links to the beginning of track
%30 and track 30 terminates without another link, then one of the entries
%of track_idxs will be [3,12,30], and the track_ids for 12 and 30 will be
%changed to 0
track_idxs = cell(ntracks,1);
track_ids = ones(ntracks,1);
for i=1:ntracks
    track_ids(i) = 1;
    track_idxs{i,1} = i;
end

for i = 1:ntracks
    idx = i;
    if ismember(idx,gaps(:,1)) && track_ids(idx)
        idxlist = idx;
        criterion = true;
        while criterion == true
            idxNext = gaps(gaps(:,1) == idx,2);
            track_ids(idxNext) = 0;
            idxlist = [idxlist,idxNext];
            if ismember(idxNext,gaps(:,1))
                idx = idxNext;
            else
                criterion = false;
            end
        end
        track_idxs{i,1} = idxlist;
    end
end

%Revise the track_idxs to get rid of duplicates using track_ids and reset
%track_ids for merging and splitting
track_idxs2 = cell(sum(track_ids),1);
idx = 1;
for i=1:size(track_idxs,1)
    if track_ids(i) == 1
        track_idxs2{idx,1} = track_idxs{i,1};
        idx = idx + 1;
    end
end
track_idxs = track_idxs2;
track_ids = ones(size(track_idxs,1),1);

%Make gap_tracks using track_idxs; this makes groups of tracks that are
%gap-closed into single tracks with discontinuities
gap_tracks = cell(length(track_ids),1);
for i=1:size(gap_tracks,1)
    gap_tracks{i} = nan(ntimepoints,1);
    for j=1:length(track_idxs{i})
        sidx = track_idxs{i}(j);
        times = simple_tracks{sidx}(:,tCol);
        idxs = simple_tracks{sidx}(:,idxCol);
        gap_tracks{i}(times) = idxs;
    end
end

%Go through merges and splits and replace track indices of simple_tracks
%with corresponding track indices for gap_tracks
splits2 = splits;
merges2 = merges;
for i = 1:size(splits2,1)
    for j = 1:size(splits2,2)-1
        sidx = splits2(i,j);
        gidx = find(cellfun(@(x) ismember(sidx,x), track_idxs));
        splits2(i,j) = gidx;
    end
end

for i = 1:size(merges2,1)
    for j = 1:size(merges2,2)-1
        sidx = merges2(i,j);
        gidx = find(cellfun(@(x) ismember(sidx,x), track_idxs));
        merges2(i,j) = gidx;
    end
end

%% Handle merging events
%Merging events are assumed to be artifacts of the detection, segmentation,
%or tracking process and should be handled before using splitting events to
%create compound tracks

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
[merge_links, bad_splits] = tempmergelinker(events, merge_cutoff);
splits2 = splits2(~ismember(splits2(:,1:2), bad_splits, 'rows'), :);
merge_links = merge_links(merge_links(:,2) ~= 0,:);

%Close gaps using merge_links; first make new track_idxs
n_gap_tracks = size(gap_tracks,1);
track_idxs = cell(n_gap_tracks,1);
track_ids = ones(n_gap_tracks,1);
for i=1:n_gap_tracks
    track_ids(i) = 1;
    track_idxs{i,1} = i;
end

%Group tracks that have been linked
for i = 1:n_gap_tracks
    idx = i;
    if ismember(idx, merge_links(:,1)) && track_ids(idx)
        idxlist = idx;
        criterion = true;
        while criterion == true
            idxNext = merge_links(merge_links(:,1) == idx,2);
            track_ids(idxNext) = 0;
            idxlist = [idxlist,idxNext];
            if ismember(idxNext, merge_links(:,1))
                idx = idxNext;
            else
                criterion = false;
            end
        end
        track_idxs{i,1} = idxlist;
    end
end

%Revise the track_idxs to get rid of duplicates using track_ids and reset
%track_ids for merging and splitting
track_idxs2 = cell(sum(track_ids),1);
idx = 1;
for i=1:size(track_idxs,1)
    if track_ids(i) == 1
        track_idxs2{idx,1} = track_idxs{i,1};
        idx = idx + 1;
    end
end
track_idxs = track_idxs2;
track_ids = ones(size(track_idxs,1),1);

%Make new_tracks using track_idxs; this makes groups of tracks that are
%gap-closed into single tracks with discontinuities
new_tracks = cell(length(track_ids),1);
for i=1:size(new_tracks,1)
    new_tracks{i} = nan(ntimepoints,1);
    for j=1:length(track_idxs{i})
        sidx = track_idxs{i}(j);
        times = find(~isnan(gap_tracks{sidx}));
        idxs = gap_tracks{sidx}(times);
        new_tracks{i}(times) = idxs;
    end
end

%Add the points from the mid track (merged to and split from) between the
%two gap-closed tracks
for k = 1:size(merge_links,1)
    for j = 1:size(track_idxs,1)
        if all(ismember(merge_links(k,1:2),track_idxs{j}))
            midtrack = merge_links(k,3);
            midtimes = merge_links(k,4):merge_links(k,5);
            midpoints = gap_tracks{midtrack}(midtimes);
            new_tracks{j}(midtimes) = midpoints;
        end
    end
end

%Go through splits2 and replace track indices of gap_tracks with 
%corresponding track indices for new_tracks
splits3 = splits2;
for i = 1:size(splits3,1)
    for j = 1:size(splits3,2)-1
        sidx = splits3(i,j);
        nidx = find(cellfun(@(x) ismember(sidx,x), track_idxs));
        splits3(i,j) = nidx;
    end
end

%% Build compound tracks using new_tracks and splits3
%make a tree object to hold compound tracks, keeping track of content and
%lineage
trees = build_tree(new_tracks, splits3);
tracks = trees;

nctracks = numel(getchildren(trees,1));
disp(strcat("Number of compound tracks = ", num2str(nctracks)))

end