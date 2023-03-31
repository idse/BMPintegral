function s =...
    split_div_cost(xy_dist,jstartA,id_A,t_dist)
%define cost for the start of track j to split from point k on track i
%ik = point k on track i
%jstart = starting point on track j
%id = daughter cell after ik if it looks like ik is dividing

diff_a = 2*abs(id_A - jstartA)./(id_A + jstartA);
s = (xy_dist + t_dist^2)*(1 + diff_a);

end