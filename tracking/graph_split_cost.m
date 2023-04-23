function s = graph_split_cost(xy_dist, ai, ai_next, aj)
%for a track midpoint i, calculate the linking cost to all track starts j
%based on position, time, area
%%% Inputs
%xyi is a 1x2 vector with the xy position of candidate cells in track i
%ti is the (scalar) frame of the midpoint in track i
%ai is the area of the last cell in track i
%xyj is an nx2 array with the xy position of all track starts within the
%cutoffs for time and distance from the midpoint of track i
%tj is an nx1 vector of midpoint times
%aj is an nx1 vector of areas of the track midpoints
%aj_prev is an nx1 vector of areas of track points immediately preceding
%these track midpoints
%%% Output
%m is an nx1 vector of costs

%denominator can be changed if needed; as it is, we are looking at the
%difference in area between ai and the sum of aj and ai_next, relative to
%aj
s = xy_dist.*(1 + abs(aj + ai_next - ai)./ai);

end