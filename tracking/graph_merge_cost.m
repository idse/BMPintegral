function m = graph_merge_cost(xyi, ai, xyj, aj, aj_prev)
%for a track end i, calculate the linking cost to all track midpoints j
%based on position, time, area
%%% Inputs
%xyi is a 1x2 vector with the xy position of candidate cells in track i
%ti is the (scalar) frame of the last point in track i
%ai is the area of the last cell in track i
%xyj is an nx2 array with the xy position of all track starts within the
%cutoffs for time and distance from the end of track i
%tj is an nx1 vector of midpoint times
%aj is an nx1 vector of areas of the track midpoints
%aj_prev is an nx1 vector of areas of track points immediately preceding
%these track midpoints
%%% Output
%m is an nx1 vector of costs

D = pdist2(xyj, xyi, 'squaredeuclidean');
%denominator can be changed if needed; as it is, we are looking at the
%difference in area between aj and the sum of ai and aj_prev, relative to
%aj
m = D(:).*(1 + abs(ai + aj_prev - aj)./aj);

end