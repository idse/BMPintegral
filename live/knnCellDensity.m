function knnCellDensity(positions,k)

ntime = positions(1).nTime;
npositions = length(positions);

for pidx = 1:npositions
    for ti = 1:ntime
        %get locations of cell centroids in the position and time selected
        XY = positions(pidx).cellData(ti).XY;
        %find squared distance between all pairs of observations and sort 
        %each column in ascending order
        D = sort(squareform(pdist(XY,'squaredeuclidean')));
        %take the kth smallest distance to other points
        d = D(k+1,:);
        %calculate the approximate cell density around each cell
        density = k./(pi*d);
        positions(pidx).cellData(ti).density = density;
    end
end


end