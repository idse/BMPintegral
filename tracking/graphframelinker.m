function [ target_indices, source_indices, total_cost, costs] =...
    graphframelinker(sinfo, tinfo, imsize, max_distance)
%The linker constructs a cost function for assigning each point in the
%source set to a point in the target set, or to nothing. 
%%Do we also want to use track history in cost computation?%%
%Get x, y, area, intensity of each point in the source and target frames:
%adjust area cutoff to determine if small enough to be dividing based
%on image magnification
%sinfo = source info, tinfo = target info

%x,y position
source = sinfo.XY;
target = tinfo.XY;
%area
source_a = sinfo.nucArea;
target_a = tinfo.nucArea;
%intensity
source_i = sinfo.nucLevel;
target_i = tinfo.nucLevel;
%marker for cell division
source_divs = sinfo.divs;
% target_divs = tinfo.divs;
%source orientations (already converted to radians and rotated by pi/2)
source_angle = sinfo.angle;

n_source_points = size(source, 1);
n_target_points = size(target, 1);

%find pairwise distance between each source and target point
D = pdist2(source,target,'squaredeuclidean');
%take advantage of MATLAB's automatic matrix expansion when doing
%arithmetic to avoid loops (column vector + row vector = matrix)
diff_i = 2*abs(target_i' - source_i)./(target_i' + source_i);
diff_a = 2*abs(target_a' - source_a)./(target_a' + source_a);
% Div = source_divs + target_divs';

%for dividing source cells, also use the orientation of the nucleus
I = ones(n_source_points, n_target_points);
idxs = find(source_divs == 1);
for ii = 1:length(idxs)
    idx = idxs(ii);
    % Pick one source point
    current_point = source(idx, :); %xy value of current point
    % Compute vector to all target points
    diff_coords = target - current_point;
    %make a unit-length vector orthogonal to the cell's major axis
    theta = source_angle(idx);
    dirvec = [cos(theta); sin(theta)];
    %normalize the distance vectors between this cell and others
    normdirs = diff_coords./sqrt(D(idx,:))';
    %use the dot product as a measure of similarity in direction
    %this ranges from 0.5 (parallel) to 1.5 (orthogonal)
    divs = 1.5 - abs(normdirs*dirvec).^3;
    I(idx,:) = divs;
end

%for pairs where we want to ignore intensity, size differences (one or both
%dividing), set diff_i, diff_a to 0
% diff_i(Div > 0) = 0;
% diff_i(Div > 1) = 0;
% diff_a(Div > 0) = 0;
% diff_a(Div > 1) = 0;
% I(Div == 1) = 1;

%combine all the costs/weights
D = D.*I;
A1 = D.*(1 + diff_i).*(1 + diff_a);
%set costs for assignments farther than the max linking distance as Inf 
A1(D > max_distance.^2) = Inf;

%Determine the alternative linking cost to be 110% of maximum cell-cell 
%linking cost (is there a better way to do this?)
linking_costs = A1(~isinf(A1)); %extract the finite costs
if ~isempty(linking_costs)
    thresh = 1.05*max(linking_costs,[],'all'); %use max value
else
    thresh = 1; %arbitrary value, linking costs are all Inf
end

%Define the weight for disappearance of each point from the source based on
%distance from edges of frame (also not the best approach)
edge_buffer = 50; %set some buffer so costs at the edge are less dramatically low
wd = 3; %make a weight to adjust the cost for disappearance
dx = min(source(:,1),imsize(1) - source(:,1)) + edge_buffer;
dy = min(source(:,2),imsize(2) - source(:,2)) + edge_buffer;
ds = wd*((dx.*dy)./(dx + dy)).^2;
ds(ds > max_distance^2) = thresh;

%Define the cost for appearance of each point in the target based on
%distance from edges of frame
wb = 3; %make a weight to adjust the cost for appearance
dx = min(target(:,1), imsize(1) - target(:,1)) + edge_buffer;
dy = min(target(:,2), imsize(2) - target(:,2)) + edge_buffer;
bs = wb*((dx.*dy)./(dx + dy)).^2;
bs(bs > max_distance^2) = thresh;

A2 = Inf(n_source_points);
A2(eye(n_source_points) == 1) = ds;

A3 = Inf(n_target_points);
A3(eye(n_target_points) == 1) = bs;

A4 = transpose(A1);

CM = [A1, A2; A3, A4];

% Find the optimal assignment using Yi Cao's optimization algorithm
[CM_indices, total_cost] = lapjv(CM); %munkres(CM);

target_indices = CM_indices(1:n_source_points);
source_indices = CM_indices((n_source_points+1):(n_source_points+n_target_points)) - n_target_points;

target_indices(target_indices > n_target_points) = -1;
source_indices(source_indices < 1) = -1;

idxs = find(target_indices > 0);

%get the associated costs for linking source points to target points
%(should be the same size as target_indices)
idxs = [idxs',target_indices(idxs)'];
ind = sub2ind(size(A1),idxs(:,1),idxs(:,2));
costs = A1(ind);



end