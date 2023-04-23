function [ target_indices, source_indices, total_cost] = livefixedlinker(source_info,...
    target_info, imsize, max_distance)
    %The linker constructs a cost function for assigning each point in the
    %source set to a point in the target set, or to nothing. 
    %Source info and target info contain columns for x, y, area, and
    %intensity of each point in the source and target frames:
    source = source_info(:,1:2);
    target = target_info(:,1:2);
    source_a = source_info(:,3);
    target_a = target_info(:,3);

    n_source_points = size(source, 1);
    n_target_points = size(target, 1);

    A1 = NaN(n_source_points, n_target_points);

    % Build distance matrix
    for i = 1 : n_source_points
        % Pick one source point
        current_point = source(i, :); %xy value of current point
        current_a = source_a(i); %area of current point
        current_a = repmat(current_a, n_target_points, 1);

        % Compute square distance to all target points
        diff_coords = target - repmat(current_point, n_target_points, 1);
        square_dist = sum(diff_coords.^2, 2);
        %Find relative difference in area between current point and target points
        diff_a = 2*abs(target_a - current_a)./(target_a + current_a); %relative difference in area
        
        % Store them
        area_weight = 1.5; %factor by which difference in area contributes to weight
        A1(i, :) = square_dist.*(1 + area_weight*diff_a);
            
    end

    % Deal with maximal linking distance: we simply mark these links as already
    % treated, so that they can never generate a link.
    A1 ( A1 > max_distance * max_distance ) = Inf;
    
    %Determine the cost to be ~>110% of maximum cell-cell linking cost
    linking_costs = A1(~isinf(A1)); %extract the finite costs
    thresh = 1.1*max(linking_costs,[],'all'); %use max value

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