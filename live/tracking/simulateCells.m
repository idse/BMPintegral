function [P, trees, sim_points] = simulateCells(numcells, ntimepoints,...
division_freq, segment_errors, velocity, intensity_drift, imsize, buffer)
disp('Making simulated data')
tic
%Point data is held in a cell array
sim_points = cell(ntimepoints,1);
%Size of image in pixels
imwidth = imsize(1);
imheight = imsize(2);
%Convert division freq to a decimal value (input should be in percent)
%Eventually do cell division in a more intelligent/realistic way
division_freq = division_freq/100;

%Create random x and y positions for the cells both inside and around the
%image frame
%Resulting x values should fall in the range (-buffer/2, imwidth + buffer/2)
%and y values should be in the range (-buffer/2, imheight + buffer/2)
Xs = (imwidth + buffer)*rand(numcells,1) - buffer/2;
Ys = (imheight + buffer)*rand(numcells,1)  - buffer/2;
%Create evenly spaced areas for the cells
As = sqrt(transpose(linspace(50,100,numcells)));
%Create random intensities for the cells (on range [0, 255] (try [10-numcells] at
%first)
Is = 29*rand(numcells,1) + 50;
%Create identifiers for the cells
identifiers = transpose(1:numcells);
%Create an index for which track this one split from - this allows
%splitting events to be traceable back to the beginning but for
%branches to still be individually identifiable - 0 is for starting track
parents = zeros(numcells,1);

%Each entry of the cell is a ncells x 3 array where each row of the array
%gives the x, y, and area values for that cell in the frame
sim_points{1} = horzcat(Xs, Ys, As, Is, identifiers, parents);

%The cells have randomly oriented motion frame-to-frame, and the order of
%the cells is randomly permuted in the array from frame-to-frame, but
%should be recoverable because each has a unique attached identifier
for i=2:ntimepoints
    ncells = size(sim_points{i-1},1);
    %Create random offset in x and y, small change in I, no change in size
    xyoff = velocity*randn(ncells,2);
    Ioff = randn(ncells,1)*intensity_drift;
    offset = horzcat(xyoff, zeros(ncells,1), Ioff, zeros(ncells,1), zeros(ncells,1));
    %Apply the offset and randomly permute the cells
    newpoint = sim_points{i-1} + offset;
    newpoint = newpoint(randperm(ncells),:);
    %randomly select some cells to divide
    dividing = rand(ncells,1) > (1 - division_freq);
    newpoint1 = newpoint(dividing,:);
    ncells1 = size(newpoint1,1);
    %set parent IDs by moving ID column to parent column
    newpoint1(:,6) = newpoint1(:,5);
    %set new track IDs
    newpoint1(:,5) = transpose(ncells+1:ncells+ncells1);
    xyoff1 = velocity*randn(ncells1,2);
    Ioff1 = randn(ncells1,1)*intensity_drift;
    offset1 = horzcat(xyoff1, zeros(ncells1,1), Ioff1, zeros(ncells1,1), zeros(ncells1,1));
    newpoint1 = newpoint1 + offset1;
    sim_points{i} = vertcat(newpoint,newpoint1);
end

%In sample_points, discard points outside the region of interest
sample_points = cell(size(sim_points));
for i=1:ntimepoints
    sampledpoints = [];
    for j=1:size(sim_points{i},1)
        jx = sim_points{i}(j,1);
        jy = sim_points{i}(j,2);
        if (jx > 0) && (jx < imwidth) && (jy > 0) && (jy < imheight)
            sampledpoints = [sampledpoints; sim_points{i}(j,:)];
        end
    end
    %Add random detection/segmentation errors
    detected_idxs = rand(size(sampledpoints,1),1) > (segment_errors/100);
    sample_points{i} = sampledpoints(detected_idxs,:);
end

%For points that overlap by a large enough margin (distance between centers
%is less than the max of the two radii), treat as one larger point with
%intensity = I1 + I2 and area = A1 + A2 - intersect(A1, A2)
overlaps = cell(size(sample_points));

%Create a new variable, sample_points2, which has an extra column allowing
%identification of one additional index; merged points have two indices,
%and non-merged points have NaN for the second index
sample_points2 = sample_points;
for i = 1:length(sample_points2)
    sample_points2{i} = [sample_points2{i}(:,1:5), NaN(size(sample_points2{i},1),1)];
end

for i=1:length(sample_points)
    merged_points = [];
    for j = 1:size(sample_points{i},1)
        center = sample_points{i}(j,1:2);
        radius = sample_points{i}(j,3);
        %dists gives the distance between the center of point i and the
        %center of other points in the frame
        dists = sample_points{i}(:,1:2) - center;
        dists = sqrt(dists(:,1).^2 + dists(:,2).^2);
        maxrad = max(radius, sample_points{i}(:,3));
        %two points are determined to be overlapping if their centers are
        %closer than the larger of the two radii
        overlap = find(dists - maxrad < 0);
        %Ignore point j overlapping with itself
        overlap = overlap(overlap ~= j);
        overlap = [j*ones(size(overlap)), overlap];
        merged_points = [merged_points; overlap];
    end
    %Since rows will generally be duplicated in merged_points (if point i
    %is overlapped with point j, then there will be a row [i,j] as well as
    %a row [j, i], these repeats should be removed:
    merged_points = unique(sort(merged_points, 2), 'rows', 'stable');
    overlaps{i} = merged_points;
    %Replace overlapping points with merged versions, where the new center,
    %(x,y) is given by ((x1 + x2)/2,(y1 + y2)/2) and the new radius is
    %given by (r1 + r2 + sqrt((x1 - x2)^2 + (y1 - y2)^2))/2
    for k = 1:size(merged_points, 1)
        idx1 = merged_points(k, 1);
        idx2 = merged_points(k, 2);
        point1 = sample_points{i}(idx1, :);
        point2 = sample_points{i}(idx2, :);
        dist = point1(1:2) - point2(1:2);
        dist = sqrt(dist(1)^2 + dist(2)^2);
        newrad = (point1(3) + point2(3) + dist)/2;
        newcenter = (point1(1:2) + point2(1:2))/2;
        newintensity = point1(4) + point2(4);
        %The first point is replaced with a merged version
        sample_points2{i}(merged_points(k,1),:) = ...
            [newcenter, newrad, newintensity, point1(5), point2(5)];
        %The second point is discarded
        sample_points2{i}(merged_points(k,2),:) = NaN;
    end
    %To remove rows of A with first element NaN, use A(~isnan(A(:,1)),:)
    sample_points2{i} = sample_points2{i}(~isnan(sample_points2{i}(:,1)),:);
end

sim_tracks = cell(size(sim_points{end},1),1);
for i = 1:ntimepoints
    for j = 1:size(sim_points{i},1)
        point = sim_points{i}(j,:);
        track_idx = point(5);
        %add the index of the point in the track and the time index
        sim_tracks{track_idx} = [sim_tracks{track_idx};[j,i]];
    end
end

%Put stuff from sample_points into a position class
P = Position();
for k = 1:ntimepoints
    P.cellData(k).XY = sample_points2{k}(:,1:2);
    P.cellData(k).area = sample_points2{k}(:,3);
    P.cellData(k).nucLevel = sample_points2{k}(:,4)*[1,1]; %This makes two identical columns
end

%Make compound ground-truth tracks using sim points
%We are interested in the last frame, because all of the tracks should be
%present in that frame
splits = sim_points{end}(:,end-1:end);
splits = splits(splits(:,2) ~= 0, :);

trees = build_tree(sim_tracks, splits);
toc
end
