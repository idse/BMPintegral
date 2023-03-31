function [ target_indices, source_indices, total_cost] = newerhungarianlinker(P,...
    time, imsize, max_distance)
%The linker constructs a cost function for assigning each point in the
%source set to a point in the target set, or to nothing. 
%%Do we also want to use track history in cost computation?%%
%Get x, y, area, intensity of each point in the source and target frames:
%adjust area cutoff to determine if small enough to be dividing based
%on image magnification

%source info
sinfo = P.cellData(time);
%target info
tinfo = P.cellData(time+1);
%x,y position
source = sinfo.XY;
target = tinfo.XY;
%area
if isfield(sinfo,'nucArea')
    source_a = sinfo.nucArea;
    target_a = tinfo.nucArea;
else
    source_a = sinfo.area;
    target_a = tinfo.area;
end
%intensity
source_i = sinfo.nucLevel(:,P.nucChannel+1);
target_i = tinfo.nucLevel(:,P.nucChannel+1);
%marker for cell division
source_divs = sinfo.divs;
target_divs = tinfo.divs;
%     target_cat = tinfo.nucMajorAxis./tinfo.nucMinorAxis;
%source orientations
source_angle = sinfo.nucOrientation*pi/180;
%     target_angle = tinfo.nucOrientation*pi/180;

n_source_points = size(source, 1);
n_target_points = size(target, 1);

A1 = NaN(n_source_points, n_target_points);

% Build distance matrix
for i = 1 : n_source_points
    % Pick one source point
    current_point = source(i, :); %xy value of current point
    current_i = source_i(i); %intensity of current point
    current_i = repmat(current_i, n_target_points, 1);
    current_a = source_a(i);
    current_a = repmat(current_a, n_target_points, 1);

    % Compute square distance to all target points
    diff_coords = target - repmat(current_point, n_target_points, 1);
    square_dist = sum(diff_coords.^2, 2);
    %Find relative difference in area & intensity between current point
    %and all target points
    diff_i = 2*abs(target_i - current_i)./(target_i + current_i);
    diff_a = 2*abs(target_a - current_a)./(target_a + current_a);

    %if the cell is small and long, reduce the cost of linking to
    %other small bright cells
%         if source_divs(i)
    %make a unit-length vector orthogonal to the cell's major axis
    theta = source_angle(i) + pi/2;
    dirvec = [cos(theta); sin(theta)];
    %normalize the distance vectors between this cell and others
%     normdirs = normalize(diff_coords,2,'norm');
%     normdirs = diff_coords./sqrt(sum(diff_coords.^2,2));
    normdirs = diff_coords./sqrt(square_dist);
    %use the dot product as a measure of similarity in direction
    %this ranges from 0.5 (parallel) to 1.5 (orthogonal)
    divs = 1.5 - abs(normdirs*dirvec).^3;
    %make the cost infinite for cells above the size cutoff
%             divs(target_a > 500) = Inf;
%             %adjust cost based on similarity in size, intensity?
%         else
%             divs = ones(n_target_points,1);
%         end
    %neither source nor target nucleus is dividing
    case1 = find(source_divs(i)+target_divs == 0);
    %only one of the source and target nuclei is dividing
    case2 = find(source_divs(i)+target_divs == 1);
    %both the source and target nuclei are dividing
    case3 = find(source_divs(i)+target_divs == 2);
    % Store them
    %A1(i, :) = square_dist.*(1 + diff_a).*(1 + diff_i);
%         A1(i, :) = square_dist.*(1 + diff_i).*divs;
    cost1 = square_dist.*(1 + diff_i).*(1 + diff_a);
    cost2 = square_dist;
    cost3 = square_dist.*divs;

    A1(i,case1) = cost1(case1);
    A1(i,case2) = cost2(case2);
    A1(i,case3) = cost3(case3);
    A1(i, square_dist > max_distance^2) = Inf;

end

% Deal with maximal linking distance: we simply mark these links as already
% treated, so that they can never generate a link.
%     A1(square_dist > max_distance^2) = Inf;

%Determine the cost to be ~>110% of maximum cell-cell linking cost
linking_costs = A1(~isinf(A1)); %extract the finite costs
if ~isempty(linking_costs)
    thresh = 1.1*max(linking_costs,[],'all'); %use max value
else
    thresh = 1; %arbitrary value, linking costs are all Inf
end

%Define the weight for disappearance of each point from the source based on
%distance from edges of frame
edge_buffer = 50; %set some buffer so costs at the edge are less dramatically low
wd = 3; %make a weight to adjust the cost for disappearance
ds = Inf(n_source_points,1);
for i=1:n_source_points
    dx = min(source(i,1), imsize(1) - source(i,1)) + edge_buffer;
    dy = min(source(i,2), imsize(2) - source(i,2)) + edge_buffer;
    tempcost = wd*(dx*dy/(dx + dy))^2;
    if tempcost > max_distance^2
        ds(i) = thresh;
    else
        ds(i) = tempcost;
    end
end

%Define the cost for appearance of each point in the target based on
%distance from edges of frame
wb = 3; %make a weight to adjust the cost for appearance
bs = Inf(n_target_points,1);
for i=1:n_target_points
    dx = min(target(i,1), imsize(1) - target(i,1)) + edge_buffer;
    dy = min(target(i,1), imsize(2) - target(i,2)) + edge_buffer;
    tempcost = wb*(dx*dy/(dx + dy))^2;
    if tempcost > max_distance^2
        bs(i) = thresh;
    else
        bs(i) = tempcost;
    end
end

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

end