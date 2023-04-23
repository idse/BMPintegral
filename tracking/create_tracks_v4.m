function G = create_tracks_v4(this, opts)
%%%% Inputs %%%%
% this                      Position object with cell data
% opts                      Structure array with following fields:
%     max_linking_distance      maximum distance for frame-frame cell 
%                               movement
%     max_gap_distance          maximum cell movement distance when 
%                               doing gap closing (distance moved over 
%                               multiple frames)
%     max_gap_closing           maximum # of frames for gap-closing, 
%                               merging, and splitting
%     imsize                    size of input image (get this from 
%                               meta instead?)
%%%% Outputs %%%%
% G                         digraph with one node for each cell in each
%                           frame and one edge for each link between cells


%parse inputs
max_linking_distance = opts.max_linking_distance;
max_gap_distance = opts.max_gap_distance;
max_gap_time = opts.max_gap_time;
imsize = opts.imsize;

if isfield(opts, 'tmax') && ~isempty(opts.tmax)
    tmax = opts.tmax;
else
    tmax = this.nTime;
end
if ~isfield(opts,'dejitter')
    opts.dejitter = false;
end
if ~isfield(opts,'divlabel')
    divlabel = 3;
else
    divlabel = opts.divlabel;
end

%replace data in empty frames
emptyframes = find(this.ncells == 0);
for ii = 1:length(emptyframes)
    frame = emptyframes(ii);
    %assumes these frames are never consecutive and that the first frame is
    %always non-empty (why would you keep imaging if the first frame is
    %bad?)
    this.cellData(frame) = this.cellData(frame-1);
end

%optionally dejitter if there is a lot of shifting of whole frames
if ~opts.dejitter
    points = transpose({this.cellData.XY});
else
    disp('Dejittering')
    [points, ~] = dejitter(this, opts.dataDir);
end

%get data on cell positions, areas, intensities, labels, orientations?
areas = transpose({this.cellData.nucArea});
%get the angle of the major axis of each cell wrt the right horizontal
%axis (i.e. positive x-axis) in radians
angles = cell(size(areas));
for ti = 1:tmax
    angles{ti} = this.cellData(ti).nucOrientation*pi/180 + pi/2;
end
%get nuclear intensity for each cell in each frame
intensities = transpose({this.cellData.nucLevel});
intensities = cellfun(@(x) x(:,this.nucChannel+1), intensities,...
    'UniformOutput', false);
%identify dividing cells
for ti = 1:tmax
    %use object classification from Ilastik to identify dividing cells
    if isfield(this.cellData,'labels')
        labels = this.cellData(ti).labels;
        this.cellData(ti).divs = labels == divlabel;
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

%initialize sparse graph with one node for each cell in each frame
ntotal = sum(this.ncells);
ntime = this.nTime;
cidxs = zeros(ntotal,2);
idx = 1;
for ti = 1:ntime
    n = this.ncells(ti);
    cidxs(idx:idx+n-1, :) = [ti*ones(n,1), (1:n)'];
    idx = idx + n;
end
%initialize a directed graph to store tracking info (this could potentially
%be done after the frame-frame linking step instead)
G = digraph(sparse(ntotal, ntotal));
%give each node in the graph an index for the timepoint and for the index
%of the cell at that timepoint
G.Nodes.frame = cidxs(:,1);
G.Nodes.cellidx = cidxs(:,2);

%%%% Frame-frame linking %%%%
disp('Doing frame-frame linking step')
tic
%Variables to hold data during linking
costs = cell(size(points));
S = cell(size(points));
T = cell(size(points));
nlinks = 0;
%for huge datasets (~10^4 cells/frame for ~500 frames) parallel for loop 
%causes the computer to run out of memory; can I write this in a smarter
%way so that doesn't happen? try not passing full positions object to the
%framelinker function each time?
parfor ti=1:tmax - 1
    sinfo = struct; %source info
    tinfo = struct; %target info
    sinfo.nucArea = areas{ti}; tinfo.nucArea = areas{ti+1};
    sinfo.XY = points{ti}; tinfo.XY = points{ti+1}; sinfo.angle = angles{ti};
    sinfo.nucLevel = intensities{ti}; tinfo.nucLevel = intensities{ti+1};
    sinfo.divs = Divs{ti}; tinfo.divs = Divs{ti+1};
    [target_idxs, ~, ~, c] = ...
        graphframelinker(sinfo, tinfo, imsize, max_linking_distance);
    costs{ti} = c(:);
    
    sidxs = find(cidxs(:,1) == ti);
    tidxs = find(cidxs(:,1) == ti + 1);
    sidxs = sidxs(target_idxs > 0);
    tidxs = tidxs(target_idxs(target_idxs > 0));
    
    nlinks = nlinks + length(sidxs);
    S{ti} = sidxs(:);
    T{ti} = tidxs(:);
    
    if mod(ti,round((tmax-1)/20)) == 0
        fprintf('.')
    end
end
fprintf('\n')

disp('adding edges to the graph')
G = addedge(G, cell2mat(S), cell2mat(T), ones(nlinks,1));
G.Edges.Cost = cell2mat(costs);
toc

%splitting, merging, gap-closing step -> build cost matrix, perform
%optimization, add new links to the graph
[assignments, costs, ~, ~] = graphtracklinker(G, max_gap_time,...
    max_gap_distance, points, areas, Divs, angles);
%put the new linking information in a convenient format for edges to be
%added to the graph
NewEdges = table(assignments, ones(size(assignments,1),1), costs, ...
    'VariableNames',{'EndNodes','Weight','Cost'});
%add new edges
G = addedge(G,NewEdges);

%resolve merging events
[newedges,badedges,nodes,costs] = newgraphmergelinker(G,points,areas,...
    Divs,max_gap_time);

obadedges = [zeros(size(nodes,1),1), nodes];
for ii = 1:length(nodes)
    obadedges(ii,1) = assignments(assignments(:,2) == nodes(ii),1);
end
badedges = [badedges; obadedges];

NewEdges = table(newedges, ones(size(newedges,1),1), costs, ...
    'VariableNames',{'EndNodes','Weight','Cost'});
if ~isempty(NewEdges)
    G = addedge(G,NewEdges);
end
if ~isempty(badedges)
    G = rmedge(G,badedges(:,1),badedges(:,2));
end


end