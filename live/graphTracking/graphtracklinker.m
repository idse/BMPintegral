function [assignments, costs, CM_size, comp_time] = ...
    graphtracklinker(G, max_gap_time, max_distance,...
    points, areas, Divs, angles)

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

%Initialize variables:
Din = indegree(G); %number of input edges to each node
Dout = outdegree(G); %number of output edges from each node
starts = find(Din == 0); %track starts have 0 in edges
ends = find(Dout == 0); %track ends have 0 out edges
midpoints = Din == 1 & Dout == 1;
ntracks = length(starts);

%These are faster/easier to index into than cell indices and correspond
%more naturally to the nodes in the graph
Frames = G.Nodes.frame;
allpoints = cell2mat(points);
allA = cell2mat(areas);
allDiv = cell2mat(Divs);
allTheta = cell2mat(angles);

%find average frame-frame displacement across all links (do this by random
%sampling instead of exhaustively)
nedges = numedges(G);
edgeids = randperm(nedges,min(nedges,1000));
avg_disp = 0;
parfor ii = 1:length(edgeids)
    nodes = G.Edges.EndNodes(edgeids(ii),:);
    xy1 = allpoints(nodes(1),:);
    xy2 = allpoints(nodes(2),:);
    avg_disp = avg_disp + sqrt(sum((xy1 - xy2).^2));
end
avg_disp = avg_disp/length(edgeids);
fprintf('average frame-frame cell displacement = %g\n', avg_disp)

%% Build cost matrices A1, A2, and B1
%Fill in the cost matrix A1 for gap closing; A(i,j) = gij = cost to close a
%gap between the end of track i and beginning of track j
disp("Building A1")
S = [];
disp(strcat("Size of A1 = ", num2str(ntracks), "x", num2str(ntracks)))
tic
parfor ii = 1:ntracks
    indi = ends(ii);
    ti = Frames(indi);
    xyi = allpoints(indi,:);
    ai = allA(indi);
    div = allDiv(indi);
    theta = allTheta(indi);
    %find indices of start tracks that are within the time and distance
    %cutoffs for gap closing
    startidxs = find(Frames(starts) > ti + 1 & Frames(starts) < ti + max_gap_time);
    if ~isempty(startidxs)
        %skip this step for speed?
        D = pdist2(xyi, allpoints(starts(startidxs),:), 'squaredeuclidean');
        startidxs = startidxs(D < max_distance^2);
    end
    %calculate the costs for linking
    if ~isempty(startidxs)
        g = graph_gap_cost(xyi, ti, ai, allpoints(starts(startidxs),:), Frames(starts(startidxs)), allA(starts(startidxs)), div, theta);
        Sijv = [ii*ones(size(g)), startidxs, g];
    else
        Sijv = [];
    end
    S = [S; Sijv];
    
    if mod(ii,round(ntracks/20)) == 0
        fprintf('.')
    end
end
fprintf('\n')
disp("Making sparse matrix")
if ~isempty(S)
    A1 = sparse(S(:,1), S(:,2), S(:,3), ntracks, ntracks);
else
    A1 = sparse(ntracks, ntracks);
end
toc

%Fill in the cost matrix A2 for merging; A2(i,j) = mij = cost to merge the
%end of track i to the middle of track j
S = [];
disp("Building A2")
tic
%For each track end (i) in simple_tracks, look at each point (k) in other
%tracks (j) that are of length 2 or greater and for which the end of i is
%temporally between the start and end of j
parfor ii = 1:ntracks
    iendidx = ends(ii); %index in allpoints
    iend = Frames(iendidx); %frame
    iendxy = allpoints(iendidx,:);
    iendA = allA(iendidx);
    %find midpoints in the timeframe of interest
    mididxs = find(midpoints & Frames == iend + 1);
    if ~isempty(mididxs)
        D = pdist2(iendxy, allpoints(mididxs,:), 'squaredeuclidean');
        mididxs = mididxs(D < max_distance^2);
    end
    
    if ~isempty(mididxs)
        aj = allA(mididxs);
        xyj = allpoints(mididxs,:);
        %find predecessors to get their areas for input; this may be faster
        %if I use find and indexing stuff instead of built-in predecessors
        %function so I can find them all at once, since I know that every
        %midpoint has exactly 1 predecessor
        aj_prev = zeros(size(aj));
        for jj = 1:length(aj_prev)
            preID = predecessors(G, mididxs(jj))
            aj_prev(jj) = allA(preID);
        end
        m = graph_merge_cost(iendxy, iendA, xyj, aj, aj_prev);
        Sijv = [ii*ones(size(mididxs)), mididxs, m];
    else
        Sijv = [];
    end
    S = [S; Sijv];
    
    if mod(ii,round(ntracks/20)) == 0
        fprintf('.')
    end
end
fprintf('\n')

disp("Making sparse matrix")
%Find and remove columns of A2 that are uniformly zero, keeping track of
%the indices of kept columns
merge_points = unique(S(:,2)); %returns column indices without repeats in ascending order
nm = length(merge_points);

for ci = 1:nm
    %this should not change size of S, only replace indices
    %ci <= keptcols(ci) always, and both are strictly increasing, so
    %nothing gets accidentally changed more than once
    S(S(:,2) == merge_points(ci),2) = ci;
end

A2 = sparse(S(:,1),S(:,2),S(:,3), ntracks, nm);
disp(strcat("Size of A2 = ", num2str(ntracks), "x", num2str(nm)))
toc

%Fill in the cost matrix B1 for splitting; B1(i,j) = sij = cost to split the
%beginning of track j from the middle of track i
S = [];
disp("Building B1")
tic
%For each midpoint (i), look at all track starts (j) and assign a cost for
%that track start to split from the midpoint
mididxs = find(midpoints);
midlength = length(mididxs);
parfor ii = 1:length(mididxs)
    mididx = mididxs(ii);
    midt = Frames(mididx);
    midxy = allpoints(mididx, :);
    mida = allA(mididx);
    div = allDiv(mididx);
    ai_next = allA(successors(G, mididx));
    
    if div
        %if midpoint is marked as dividing, get position, size of daughter
        %cell in subsequent frames
        id_t = (midt+1:midt + max_gap_time)';
        id_xy = zeros(max_gap_time,2);
        id_A = zeros(max_gap_time,1);
        count = 1;
        criterion = true;
        nodeID = mididx;
        while criterion && count < max_gap_time
            nodeID = successors(G, nodeID);
            if isempty(nodeID)
                criterion = false;
            else
                id_xy(count,:) = allpoints(nodeID,:);
                id_A(count) = allA(nodeID);
            end
            count = count + 1;
        end
        id_t = id_t(1:count-1);
        id_xy = id_xy(1:count-1,:);
        id_A = id_A(1:count-1);
        %expected position for a daughter cell is based on position of
        %first daughter cell
        ex_xy = 2*midxy - id_xy;
        maxt = id_t(end);
        %find start tracks in the appropriate time window
        sidxs = find(Frames(starts) > midt & Frames(starts) == midt + max_gap_time);
        startidxs = starts(sidxs);
        if ~isempty(startidxs)
            jxy = allpoints(startidxs,:);
            jt = Frames(startidxs);
            ja = allA(startidxs);
            %find the appropriate position and area against which to compare
            %each track start
            ixy = zeros(size(jxy));
            ia = zeros(size(ja));
            for ti = 1:length(id_t)
                times = find(jt == id_t(ti));
                if ~isempty(times)
                    ixy(times,:) = id_xy(ti,:);
                    ia(times) = id_A(ti);
                end
            end
            xy_dist = sum((jxy - ixy).^2,2);
            sij = graph_split_div_cost(xy_dist, ja, ia);
            idxs = xy_dist < max_distance;
            sij = sij(idxs);
            startidxs = startidxs(idxs);
            sidxs = sidxs(idxs);
            if ~isempty(sij)
                disp('div link cost added')
                Sijv = [mididx*ones(size(startidxs)), sidxs, sij];
            else
                Sijv = [];
            end
        else
            Sijv = [];
        end
    else
        %if not dividing
        %indices just among the starts
        sidxs = find(Frames(starts) == midt + 1);
        %indices among all points
        startidxs = starts(sidxs);
        xyj = allpoints(startidxs,:);
        if ~isempty(startidxs)
            xy_dist = pdist2(midxy, xyj, 'squaredeuclidean')';
            idxs = xy_dist < max_distance^2;
            startidxs = startidxs(idxs);
            sidxs = sidxs(idxs);
        end
        %if there are start points in the right time frame within the
        %cutoff distance of this point, calculate the costs
        if ~isempty(startidxs)
            sij = graph_split_cost(xy_dist(idxs), mida, ai_next, allA(startidxs));
            Sijv = [mididx*ones(size(sidxs)), sidxs, sij];
        else
            Sijv = [];
        end
    end
    S = [S; Sijv];
    
    if mod(ii,round(length(midlength)/20)) == 0
        fprintf('.')
    end
end
fprintf('\n')

disp("Making sparse matrix")
%Find and remove rows of B1 that are uniformly zero, keeping track of the
%indices of kept rows
split_points = unique(S(:,1)); %returns unique row indices in sorted order
ns = length(split_points);
for ri = 1:ns
    S(S(:,1) == split_points(ri),1) = ri; %ignore this warning
end

B1 = sparse(S(:,1),S(:,2),S(:,3), ns, ntracks);
disp(strcat("Size of B1 = ", num2str(ns), "x", num2str(ntracks)))
toc


%% Build alternative cost matrices D2 and B4 to reject merging and splitting
%Build matrix A5 with alternative costs to reject merging events
disp('Building B4 and D2')
%make separate functions for the alternative linking costs
tic
%alternative cost matrix for splitting
altcosts = zeros(ns, 1);
for jj = 1:ns
    aj = allA(split_points(jj));
    aj_prev = allA(predecessors(G,split_points(jj)));
    altcosts(jj) = avg_disp^2*(1 + abs(aj - aj_prev)/aj);
end
B4 = spdiags(altcosts, 0, ns, ns);

%alternative cost matrix for merging
altcosts = zeros(nm, 1);
for ii = 1:nm
    ai = allA(merge_points(ii));
    ai_next = allA(successors(G, merge_points(ii)));
    altcosts(ii) = avg_disp^2*(1 + abs(ai - ai_next)/ai);
end
D2 = spdiags(altcosts, 0, nm, nm);
toc

%% Build alternative cost matrices A3 and C1 to reject gap closing
%Gap closing rejection cost is based on proximity to image boundary, but
%scaled up? Doesn't work because cells close to edge of frame are more
%likely to have disappeared for one or more frames then reappeared so you
%dont want to reject these

%Preliminary cost: assign constant value that is in ~110th percentile of 
%gap closing costs
gapcosts = sort(nonzeros(A1));
if length(gapcosts) > 1
    %make gap-rejection costs 1.1*max non-inf gap-closing cost
    bscale = 1.1*gapcosts(end);
else
    bscale = 1;
end

A3 = spdiags(bscale*ones(ntracks,1), 0, ntracks, ntracks);
C1 = spdiags(bscale*ones(ntracks,1), 0, ntracks, ntracks);

%% Define infinite and duplicate (transpose) matrices
% This is just to meet to requirements of the mathematical formalism
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
CM = [A1, A2, A3, A4;...
    B1, B2, B3, B4;...
    C1, C2, C3, C4;...
    D1, D2, D3, D4];

%find the maximum non-zero (non-inf) cost in the cost matrix and scale it
%by the height of the matrix to get a very large offset relative to the
%values in the matrix
off = size(CM,1)*max(nonzeros(CM));
%apply the offset to the non-zero elements of the cost matrix so they are
%all very negative and costs of zero are relatively huge
CM = spfun(@(x) x - off, CM);
CM_size = size(CM);
whos CM

disp('Performing optimization')
tic
[track_indices, ~] = lapjv_v2(CM);
comp_time = toc;
disp(strcat("Elapsed time is ", num2str(comp_time), " seconds"))

%% Parse the output indices to assign closing, merging, splitting events
%determine the start and end node in the graph for each assigned (not
%rejected) link

%make costs positive again so that they can be extracted from the matrix
%for a comparison of relative confidences in assignments later
CM = spfun(@(x) x + off, CM);

gap_counts = 0;
merge_counts = 0;
split_counts = 0;

%assignments(jj,:) = [start_node, end_node];
assignments = NaN(ntracks + ns, 2);
costs = NaN(ntracks + ns, 1);
%gap closing and merging
for ii = 1:ntracks
    idx = track_indices(ii);
    endidx = ends(ii);
    if idx <= ntracks
        gap_counts = gap_counts + 1;
        assignments(ii, :) = [endidx, starts(idx)];
        costs(ii) = CM(ii, idx);
    elseif idx <= ntracks + nm
        merge_counts = merge_counts + 1;
        mid_point = merge_points(idx-ntracks);
        assignments(ii, :) = [endidx, mid_point];
        costs(ii) = CM(ii, idx);
    end
end
%splitting
for ii = (ntracks + 1):(ntracks + ns)
    mid_point = split_points(ii - ntracks);
    idx = track_indices(ii);
    if idx <= ntracks
        split_counts = split_counts + 1;
        assignments(ii, :) = [mid_point, starts(idx)];
        costs(ii) = CM(ii, idx);
    end
end

disp(strcat("Number of gaps closed = ", num2str(gap_counts)))
disp(strcat("Number of merges = ", num2str(merge_counts)))
disp(strcat("Number of splits = ", num2str(split_counts)))

nanidxs = isnan(costs);
assignments(nanidxs,:) = [];
costs(nanidxs) = [];

end