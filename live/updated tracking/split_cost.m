function s =...
    split_cost(xy_dist,ik_t,ik_A,jstart,jstartA, inext_A)
%define cost for the start of track j to split from point k on track i
%ik = point k on track i
%jstart = starting point on track j
%id = daughter cell after ik if it looks like ik is dividing

t_dist = jstart - ik_t;
pij = ik_A/(jstartA + inext_A);
if pij >= 1
    s = (xy_dist + t_dist^2)*pij;
elseif pij < 1 && pij > 0
    s = (xy_dist + t_dist^2)/pij^2;
elseif pij <= 0
    disp('Something is wrong')
end

end