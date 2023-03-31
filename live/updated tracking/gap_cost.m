function g = gap_cost(xj,xi,yj,yi,tj,ti,aj,ai,div,theta)
%calculate the cost to link the end of track i to the beginning of track j
%given the coordinates, time point, and area of the nuclei at the end of
%track i and the beginning of track j
%if the nucleus at the end of track i appears to be dividing, define the
%cost in a way that reflects the behavior we expect from a dividing nucleus

if div
    xyvec = [xj - xi; yj - yi];
    xyvec = xyvec/norm(xyvec);
    dirvec = [cos(theta); sin(theta)];
    div_weight = 1.5 - abs(xyvec'*dirvec)^3;
    g = div_weight*((xj - xi)^2 + (yj - yi)^2 + (tj - ti)^2);
else
    g = ((xj - xi)^2 + (yj - yi)^2 + (tj - ti)^2)*(1 + 2*abs(ai - aj)/(ai + aj));
end

end