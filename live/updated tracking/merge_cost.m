function m = merge_cost(xy_dist, t_dist, iendA, jk_A, jprev_A)

pij = jk_A/(iendA + jprev_A);
if pij > 1
    m = (xy_dist + t_dist^2)*pij;
elseif pij <= 1 && pij > 0
    m = (xy_dist + t_dist^2)/pij^2;
elseif pij <= 0
    disp('Something is wrong')
end

end