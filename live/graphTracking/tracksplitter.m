function [assignments, costs, CM_size, comp_time] = ...
    tracksplitter(G, max_gap_time, max_distance,...
    points, areas, Divs)

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
midpoints = Din == 1 & Dout == 1; %midpoints have edge in and 1 edge out
ntracks = length(starts);

%These are faster/easier to index into than cell indices and correspond
%more naturally to the nodes in the graph
Frames = G.Nodes.frame;
allpoints = cell2mat(points);
allA = cell2mat(areas);
allDiv = cell2mat(Divs);

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

%% Build cost matrix A1 for splitting
%Fill in the cost matrix A1 for splitting; A1(i,j) = sij = cost to split the
%beginning of track j from the middle of track i
S = [];
disp("Building A1")
tic
%For each midpoint (i), look at all track starts (j) and assign a cost for
%that track start to split from the midpoint
mididxs = find(midpoints);
midlength = length(mididxs);
for ii = 1:length(mididxs)
    mididx = mididxs(ii);
    midt = Frames(mididx);
    midxy = allpoints(mididx, :);
%     mida = allA(mididx);
    div = allDiv(mididx);
%     ai_next = allA(successors(G, mididx));
    
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
        sidxs = find(Frames(starts) > midt & Frames(starts) <= midt + max_gap_time);
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
%                     ixy(times,:) = ex_xy(ti,:);%id_xy(ti,:);
                    ixy(times,1) = ex_xy(ti,1);
                    ixy(times,2) = ex_xy(ti,2);
                    ia(times) = id_A(ti);
                end
            end
            xy_dist = sum((jxy - ixy).^2,2);
            sij = graph_split_div_cost(xy_dist, ja, ia);
            idxs = xy_dist < max_distance;
            sij = sij(idxs);
            sidxs = sidxs(idxs);
            startidxs = startidxs(idxs);
            if ~isempty(sij)
                Sijv = [mididx*ones(size(startidxs)), sidxs, sij];
            else
                Sijv = [];
            end
        else
            Sijv = [];
        end
    else
        %don't split from non-dividing nuclei
        Sijv = [];
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
if ~isempty(S)
    split_points = unique(S(:,1)); %returns unique row indices in sorted order
    ns = length(split_points);
    for ri = 1:ns
        S(S(:,1) == split_points(ri),1) = ri; %ignore this warning
    end

    A1 = sparse(S(:,1),S(:,2),S(:,3), ns, ntracks);
else
    %if there are no possible splitting events, just return
    %is there a better way to do this? i.e. set ns = 0, A1 = []
    assignments = [];
    costs = [];
    CM_size = [0 0];
    comp_time = toc;
    return
end
disp(strcat("Size of A1 = ", num2str(ns), "x", num2str(ntracks)))
disp(strcat("Number of nonzero elements = ", num2str(nnz(A1))))
toc


%% Build alternative cost matrix A2 to reject splitting
%Build matrix A5 with alternative costs to reject merging events
disp('Building B4 and D2')
%make separate functions for the alternative linking costs
tic
%alternative cost matrix for splitting
altcosts = zeros(ns, 1);
for jj = 1:ns
    aj = allA(split_points(jj));
    aj_prev = allA(predecessors(G,split_points(jj)));
    altcosts(jj) = 3*avg_disp^2*(1 + abs(aj - aj_prev)/aj);
end
A2 = spdiags(altcosts, 0, ns, ns);
toc

%% Build alternative cost matrix B1 for track initiation
%Gap closing rejection cost is based on proximity to image boundary, but
%scaled up? Doesn't work because cells close to edge of frame are more
%likely to have disappeared for one or more frames then reappeared so you
%dont want to reject these

%Preliminary cost: assign constant value that is in ~110th percentile of 
%splitting costs
gapcosts = sort(nonzeros(A1));
if length(gapcosts) > 1
    %make gap-rejection costs 1.1*max non-inf gap-closing cost
    bscale = 1.1*gapcosts(end);
else
    bscale = 1;
end

B1 = spdiags(bscale*ones(ntracks,1), 0, ntracks, ntracks);

%% Construct full cost matrix and perform optimization
B2 = transpose(A1);

CM = [...
    A1, A2;...
    B1, B2];

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

%% Parse the output indices to assign splitting events
%determine the start and end node in the graph for each assigned (not
%rejected) link

%make costs positive again so that they can be extracted from the matrix
%for a comparison of relative confidences in assignments later
CM = spfun(@(x) x + off, CM);

split_counts = 0;

%assignments(jj,:) = [start_node, end_node];
assignments = NaN(ns, 2);
costs = NaN(ns, 1);
%splitting
for ii = 1:ns
    mid_point = split_points(ii);
    idx = track_indices(ii);
    if idx <= ntracks
        split_counts = split_counts + 1;
        assignments(ii, :) = [mid_point, starts(idx)];
        costs(ii) = CM(ii, idx);
    end
end

disp(strcat("Number of splits = ", num2str(split_counts)))

nanidxs = isnan(costs);
assignments(nanidxs,:) = [];
costs(nanidxs) = [];

end