function [gaps, merges, splits, total_cost, CM_size, comp_time] = tracklinker8(simple_tracks,...
    max_gap_distance, max_distance, points, intensities)

%We will build a 4x4 block cost matrix:
%  | A1   A2   A3   Inf |
%  | B1   Inf  Inf  B4  |
%  | C1   Inf  C3   C4  |
%  | Inf  D2   D3   Inf |

%Where blocks correspond to:
%A1 - gap closing
%A2 - merging
%A3 - terminating
%B1 - splitting
%B4 - no splitting
%C1 - initiating
%C3 = transpose(A1) - gap closing
%C4 = transpose(B1) - splitting
%D2 - no merging
%D3 = transpose(A2) - merging

%And then do the optimization and parse the outputs to determine where
%there is gap closing, merging, splitting, terminating, and initiating

%suppress warning about temporary variables in parfor loop:
warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary')

disp("Using sparse matrix formulation")

%Initialize variables:
ntracks = size(simple_tracks, 1);
ntime = length(points);
track_lengths = cellfun(@(x) size(x,1), simple_tracks);
tracklentotal = sum(track_lengths);
disp(strcat("Cumulative track length = ", num2str(tracklentotal)))

tCol = 2;
idxCol = 1;

%Store the index and time point for the beginning and end of each track
track_starts = zeros(ntracks,2);
track_ends = zeros(ntracks,2);
for k=1:ntracks
    track_starts(k,:) = simple_tracks{k}(1,:);
    track_ends(k,:) = simple_tracks{k}(end,:);
end
%summed_lengths gives the sum of the lengths of all the tracks (should be
%equal to the total number of points across all frames)
summed_lengths = cumsum(track_lengths);

%Find average frame-frame displacements for each track segment; these will
%be used in constructing alternative costs to reject merging and splitting
avg_disps = zeros(length(simple_tracks),1);
for k = 1:length(simple_tracks)
    total_disp = 0;
    %If a track segment only has one point, the average displacement is 0
    if size(simple_tracks{k},1) > 1
        for i=2:size(simple_tracks{k},1)
            idx = simple_tracks{k}(i,idxCol);
            idxlast = simple_tracks{k}(i-1,idxCol);
            time = simple_tracks{k}(i,tCol);
            total_disp = total_disp + sqrt(sum((points{time-1}(idxlast,1:2) - points{time}(idx,1:2)).^2));
        end
        avg_disp = total_disp/(length(simple_tracks{k})-1);
        avg_disps(k) = avg_disp;
    end
end

%Will need a good way to keep track of which tracks/points are actually in
%the cost matrices for merging and splitting
merge_points = zeros(sum(track_lengths), 2);
merge_points(1:track_lengths(1),1:2) = ...
        [(1:track_lengths(1))', ones(track_lengths(1),1)];
for k = 2:ntracks
    merge_points(summed_lengths(k-1)+1:summed_lengths(k), 1:2) = ...
        [(1:track_lengths(k))', k*ones(track_lengths(k), 1)];
end

split_points = zeros(sum(track_lengths), 2);
split_points(1:track_lengths(1),1:2) = ...
        [(1:track_lengths(1))', ones(track_lengths(1),1)];
for k = 2:ntracks
    split_points(summed_lengths(k-1)+1:summed_lengths(k), 1:2) = ...
        [(1:track_lengths(k))', k*ones(track_lengths(k), 1)];
end

%% Build cost matrices A1, A2, and B1
%Fill in the cost matrix A1 for gap closing; A(i,j) = gij = cost to close a
%gap between the end of track i and beginning of track j
disp("Building A1")
S = [];
disp(strcat("Size of A1 = ", num2str(ntracks), "x", num2str(ntracks)))
tic
parfor i=1:ntracks
    Sijv = zeros(ntracks,3);
    ti = simple_tracks{i}(end,2);
    indi = simple_tracks{i}(end,1);
    xi = points{ti}(indi,1);
    yi = points{ti}(indi,2);
    for j=1:ntracks
        tj = simple_tracks{j}(1,2);
        indj = simple_tracks{j}(1,1);
        xj = points{tj}(indj,1);
        yj = points{tj}(indj,2);
        xydist = (xi - xj)^2 + (yi - yj)^2;
        tdist = tj - ti;
        %If the frame-to-frame distance is greater than the max gap closing
        %or if the x-y distance is greater than the max linking distance,
        %set the cost = Inf; otherwise cost = delx^2 + dely^2 + delt^2
        if (tdist > 2) && (tdist < max_gap_distance) && (xydist <= max_distance^2) && (i ~= j)
            gij = (xj - xi)^2 + (yj - yi)^2 + (tj - ti)^2;
            Sijv(j,:) = [i, j, gij];
        end
    end
    S = [S; Sijv(Sijv(:,1) ~= 0,:)];
    
    if mod(i,round(ntracks/20)) == 0
        disp(".")
    end
end
disp("Making sparse matrix")
A1 = sparse(S(:,1), S(:,2), S(:,3), ntracks, ntracks);
toc

%Fill in the cost matrix A2 for merging; A2(i,j) = mij = cost to merge the
%end of track i to the middle of track j
S = [];
disp("Building A2")
% disp(strcat("Size of A2 = ", num2str(ntracks), "x", num2str(tracklentotal)))
tic
%For each track end (i) in simple_tracks, look at each point (k) in other
%tracks (j) that are of length 2 or greater and for which the end of i is
%temporally between the start and end of j
parfor i = 1:ntracks
    Sijv = zeros(tracklentotal,3);
    iend = track_ends(i,tCol);
    iendidx = track_ends(i,idxCol);
    iendx = points{iend}(iendidx,1);
    iendy = points{iend}(iendidx,2);
    iendInt = intensities{iend}(iendidx);
    for j = 1:ntracks
        jstart = track_starts(j, tCol);
        jend = track_ends(j, tCol);
        %If track j is of length 2 or greater and iend is in [jstart,jend)
        %iend<jend should preclude cases where i = j (can't merge the same
        %track to itself)
        if (track_lengths(j) > 1) && (iend >= jstart) && (iend < jend)
            jtimes = simple_tracks{j}(:,tCol);
            %Only loop over relevant points - this finds the indices of
            %points with time > end of i and within the max gap closing
            %window (and doesn't allow merging in the last frame)
            rel_points = intersect(find(jtimes > iend), find(jtimes < min([iend + max_gap_distance, ntime, jend])));
            for l = 2:length(rel_points)
                k = rel_points(l);
                %Pull out time, x, y, and intensity values for point k on
                %track j in order to calculate the linking cost
                jk_idx = simple_tracks{j}(k,idxCol);
                jk_t = jtimes(k);
                jk_x = points{jk_t}(jk_idx,1);
                jk_y = points{jk_t}(jk_idx,2);
                jk_int = intensities{jk_t}(jk_idx);
                xy_dist = (jk_x - iendx)^2 + (jk_y - iendy)^2;
                t_dist = jk_t - iend;
                if xy_dist <= max_distance^2
                    jprev_idx = simple_tracks{j}(k-t_dist,idxCol);
                    jprev_int = intensities{iend}(jprev_idx);
                    pij = jk_int/(iendInt + jprev_int);
                    if pij > 1
                        mij = (xy_dist + t_dist^2)*pij;
                    elseif pij < 1 && pij > 0
                        mij = (xy_dist + t_dist^2)/pij^2;
                    elseif pij <= 0
                        disp('Something is wrong')
                    end
                    jidx = k + sum(track_lengths(1:j-1));
                    Sijv(jidx,:) = [i, jidx, mij];
                end
            end
        end
    end
    S = [S; Sijv(Sijv(:,1) ~= 0,:)];
    
    if mod(i,round(ntracks/20)) == 0
        disp('.')
    end
end

disp("Making sparse matrix")
%Find and remove columns of A2 that are uniformly zero, keeping track of
%the indices of kept columns
keptcols = unique(S(:,2)); %returns column indices without repeats in ascending order
merge_points = merge_points(keptcols,:);
nm = size(merge_points,1);

for ci = 1:nm
    S(S(:,2) == keptcols(ci),2) = ci; %this should not change size of S, only replace indices
end

A2 = sparse(S(:,1),S(:,2),S(:,3), ntracks, nm);
disp(strcat("Size of A2 = ", num2str(ntracks), "x", num2str(nm)))
% A2 = sparse(S(:,1), S(:,2), S(:,3), ntracks, tracklentotal);
toc

%Fill in the cost matrix B1 for splitting; B1(i,j) = sij = cost to split the
%beginning of track j from the middle of track i
% B1 = cell(1, ntracks); %zeros(tracklentotal, ntracks);
S = [];
disp("Building B1")
% disp(strcat("Size of B1 = ", num2str(tracklentotal), "x", num2str(ntracks)))
tic
%For each track start (j) in simple_tracks, look at each point (k) in other
%tracks (i) that are of length 2 or greater and for which the start of j is
%temporally between the start and end of i
% B1height = size(B1,1);
parfor j = 1:ntracks
    Sijv = zeros(tracklentotal, 3);
    
    jstart = track_starts(j,tCol);
    jstartidx = track_starts(j,idxCol);
    jstartx = points{jstart}(jstartidx,1);
    jstarty = points{jstart}(jstartidx,2);
    jstartInt = intensities{jstart}(jstartidx);
    for i=1:ntracks
        istart = track_starts(i,tCol);
        iend = track_ends(i,tCol);
        %If track i is of length 2 or greater and jstart is in (istart,iend]
        %jstart>istart should preclude cases where i = j (can't split the same
        %track from itself)
        if (track_lengths(i) > 1) && (jstart <= iend) && (jstart > istart)
            itimes = simple_tracks{i}(:,tCol);
            %Only loop over relevant points - this finds the indices of
            %points with time < start of j and within the max gap closing
            %window
            rel_points = intersect(find(itimes < jstart), find(itimes > jstart - max_gap_distance));
            for l = 2:length(rel_points)
                k = rel_points(l);
                %Extract time, x, y, and intensity values for point k on
                %track i in order to calculate the linking cost
                ik_idx = simple_tracks{i}(k,idxCol);
                ik_t = itimes(k);
                ik_x = points{ik_t}(ik_idx,1);
                ik_y = points{ik_t}(ik_idx,2);
                ik_int = intensities{ik_t}(ik_idx);
                xy_dist = (ik_x - jstartx)^2 + (ik_y - jstarty)^2;
                t_dist = jstart - ik_t;
                if xy_dist <= max_distance^2
                    inext_idx = simple_tracks{i}(k + t_dist,idxCol);
                    inext_int = intensities{jstart}(inext_idx);
                    pij = ik_int/(jstartInt + inext_int);
                    if pij > 1
                        sij = (xy_dist + t_dist^2)*pij;
                    elseif pij < 1 && pij > 0
                        sij = (xy_dist + t_dist^2)/pij^2;
                    elseif pij <= 0
                        disp('Something is wrong')
                    end
                    iidx = k + sum(track_lengths(1:i-1));
                    Sijv(iidx,:) = [iidx, j, sij];
                end
            end
        end
    end
    S = [S; Sijv(Sijv(:,1) ~= 0,:)];
    
    if mod(j,round(ntracks/20)) == 0
        disp('.')
    end
end
disp("Making sparse matrix")
%Find and remove rows of B1 that are uniformly zero, keeping track of the
%indices of kept columns
keptrows = unique(S(:,1)); %returns unique row indices in sorted order
split_points = split_points(keptrows,:);
ns = size(split_points,1);

for ri = 1:ns
    S(S(:,1) == keptrows(ri),1) = ri; %ignore this warning
end

B1 = sparse(S(:,1),S(:,2),S(:,3), ns, ntracks);
disp(strcat("Size of B1 = ", num2str(ns), "x", num2str(ntracks)))
% B1 = sparse(S(:,1), S(:,2), S(:,3), tracklentotal, ntracks);
toc

% %Find points for which both the corresponding column of A2 and row of B1
% %are uniformly infinite-valued; these points can be disregarded
% infcols = find(all(A2 == 0, 1));
% infrows = find(all(B1 == 0, 2));
% keptcols = find(~ismember(1:tracklentotal, infcols));
% keptrows = find(~ismember(1:tracklentotal, infrows));
% 
% %Get rid of rows and columns of B1 and A2 that may be ignored, keeping
% %track of which ones are kept
% A2 = A2(:,keptcols);
% B1 = B1(keptrows,:);
% merge_points = merge_points(keptcols,:);
% split_points = split_points(keptrows,:);
% ns = size(split_points,1);
% nm = size(merge_points,1);

%% Build alternative cost matrices D2 and B4 to reject merging and splitting
%Build matrix A5 with alternative costs to reject merging events
%The alternative cost for some point k in track segment j is given by
%bm = (delbar_j)^2 * pij, where delbar_j is the average frame-frame
%displacement for track segment j and pij is Int_j(t)/Int_j(t-1) (if pij <
%1, then invert it)
disp('Building B4 and D2')
tic
altcosts = zeros(tracklentotal, 1);
current = 1;
unitcost = 1000; %Some arbitrary alternative cost for tracks of length 1

%This computes tons of unused costs - later make it only iterate over keep,
%using merge_points to figure out which track and point to use to find the
%cost (but takes less than a second now, so efficiency isn't a pressing
%concern)
for i = 1:length(simple_tracks)
    if mod(i,ceil(length(simple_tracks)/20)) == 0
        disp(strcat(num2str(round(i*100/length(simple_tracks))),'%'))
    end
    avg_disp = avg_disps(i);
    itrack = simple_tracks{i};
    for j = 1:size(itrack,1)
        %If track is of length 1, assign unitcost
        if size(itrack,1) == 1
            altcosts(current) = unitcost;
        %If last point in track, assign same cost as previous point
        elseif j == size(itrack,1)
            altcosts(current) = altcosts(current - 1);
        %Otherwise, base cost on average frame-frame displacement and the
        %relative change in intensity for this point from the previous one
        else
            idx_j = itrack(j,idxCol);
            t_j = itrack(j,tCol);
            idxNext = itrack(j+1,idxCol);
            I_t = intensities{t_j}(idx_j);
            I_next = intensities{t_j + 1}(idxNext);
            pj = max(I_t/I_next, I_next/I_t);
            altcosts(current) = pj*avg_disp^2;
        end
        current = current + 1;
    end
end

B4 = spdiags(altcosts, 0, size(split_points,1), size(split_points,1));
D2 = spdiags(altcosts, 0, size(merge_points,1), size(merge_points,1));

toc

%% Build alternative cost matrices A3 and C1 to reject gap closing
%Gap closing rejection cost is based on proximity to image boundary, but
%scaled up? Doesn't work because cells close to edge of frame are more
%likely to have disappeared for one or more frames then reappeared

%Preliminary cost: assign constant value that is in ~90th percentile of gap
%closing costs
% gapcosts = sort(A1(A1 ~= 0));
gapcosts = sort(nonzeros(A1));
if length(gapcosts) > 1
    bscale = gapcosts(round(0.9*length(gapcosts)));
else
    bscale = 1;
end

A3 = spdiags(bscale*ones(ntracks,1), 0, ntracks, ntracks);
C1 = spdiags(bscale*ones(ntracks,1), 0, ntracks, ntracks);

%% Define infinite and duplicate (transpose) matrices
% This is just to meet to topological requirements of the mathematical
% formalism
A4 = sparse(ntracks, ns);
B2 = sparse(ns, nm);
B3 = sparse(ns, ntracks);
C2 = sparse(ntracks, nm);
D1 = sparse(nm, ntracks);
D4 = sparse(nm, ns);
C3 = transpose(A1);
C4 = transpose(B1);
D3 = transpose(A2);


%% Construct full cost matrix and perform optimization
%CM = [A1, A2, A3; B1, A5, A2'; C1, B1', A1'];
%try storing everything in sparse matrices up until the actual optimization
%must be performed, then converting to a dense matrix and replacing zeros
%with Infs
CM = [A1, A2, A3, A4;...
    B1, B2, B3, B4;...
    C1, C2, C3, C4;...
    D1, D2, D3, D4];

off = size(CM,1)*max(nonzeros(CM));
CM = spfun(@(x) x - off, CM);
CM_size = size(CM);
whos CM

disp('Performing optimization')
tic
[track_indices, total_cost] = lapjv_v2(CM);
comp_time = toc;
disp(strcat("Elapsed time is ", num2str(comp_time), " seconds"))

%% Parse the output indices to assign closing, merging, splitting events
%Maybe add a double check that all the indices are consistent and the
%number of each event is counted the same way both ways
%gaps(i,:) = [end_track, start_track, gap_time, cost]
%merges(i,:) = [end_track, mid_track, merge_time, cost]
%splits(i,:) = [start, mid_track, split_time, cost]

%make costs positive again so that they can be extracted from the matrix
%for a comparison of relative confidences in assignments later
CM = spfun(@(x) x + off, CM);

gap_counts = 0;
merge_counts = 0;
split_counts = 0;

gaps = -1*ones(ntracks, 4);
merges = -1*ones(ntracks, 4);
splits = -1*ones(ntracks, 4);

for i = 1:ntracks
    idx = track_indices(i);
    if idx <= ntracks
        gap_counts = gap_counts + 1;
        gap_time = track_starts(idx, tCol);
        gaps(gap_counts,:) = [i, idx, gap_time, CM(i,idx)];
    elseif idx <= ntracks + nm
        merge_counts = merge_counts + 1;
        mid_track = merge_points(idx-ntracks,2);
        mid_point = merge_points(idx-ntracks,1);
        merge_time = simple_tracks{mid_track}(mid_point, tCol);
        merges(merge_counts,:) = [i, mid_track, merge_time, CM(i,idx)];
    end
end

for i = (ntracks + 1):(ntracks + ns)
    idx = track_indices(i);
    if idx <= ntracks
        split_counts = split_counts + 1;
        mid_track = split_points(i-ntracks, 2);
        mid_point = split_points(i-ntracks, 1);
        split_time = simple_tracks{mid_track}(mid_point, tCol);
        splits(split_counts,:) = [idx, mid_track, split_time, CM(i,idx)];
    end
end

disp(strcat("Number of gaps closed = ", num2str(gap_counts)))
disp(strcat("Number of merges = ", num2str(merge_counts)))
disp(strcat("Number of splits = ", num2str(split_counts)))

gaps(all(gaps == -1,2),:) = [];
merges(all(merges == -1,2),:) = [];
splits(all(splits == -1,2),:) = [];

%turn the warning back on at the end
warning('on', 'MATLAB:mir_warning_maybe_uninitialized_temporary')

end