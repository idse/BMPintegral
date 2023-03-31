function s = graph_split_div_cost(xy_dist, ja, ia)
%jxy is an nx2 array with the x-y coordinates of track starts
%ja is an nx1 vector of the areas of track start cells
%ixy is an nx2 array with the expected xy position against which to compare
%each of these track starts if it split from the midpoint of track i (based
%on the time indices of the points in jxy)
%ia is, likewise, an nx1 vector of expected track areas

s = xy_dist.*(1 + 2*abs(ia - ja)./(ia + ja));


end