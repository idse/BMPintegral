function g = graph_gap_cost(xyi, ti, ai, xyj, tj, aj, div, theta)
%for a track end i, calculate the linking cost to all track starts j based
%on position, time, area
%%% Inputs
%xyi is a 1x2 vector with the xy position of the last cell in track i
%ti is the (scalar) frame of the last point in track i
%ai is the area of the last cell in track i
%xyj is an nx2 array with the xy position of all track starts within the
%cutoffs for time and distance from the end of track i
%tj is an nx1 vector of start times for those tracks
%aj is an nx1 vector of areas of the first points of those tracks
%%% Output
%g is an nx1 vector of costs

D = pdist2(xyj, xyi, 'squaredeuclidean');

if div
    diff_coords = xyj - xyi;
    normdirs = diff_coords./sqrt(D);
    dirvec = [cos(theta); sin(theta)];
    div_weight = 1.5 - abs(normdirs*dirvec).^3;
    g = div_weight.*(D + (tj - ti).^2);
else
    g = (D + (tj - ti).^2).*(1 + abs(ai - aj)./ai);
end


end