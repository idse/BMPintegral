function colorscatter(x,y,z,tol,pointsize)

if ~exist('tol','var')
    tol = 0.02;
end

if ~exist('pointsize','var')
    pointsize = 15;
end

nnz = z(~isnan(z));
n = length(nnz);
zs = sort(nnz);
zmin = zs(max(ceil(n*tol),1));
zmax = zs(floor(n*(1-tol)));
z(z < zmin) = zmin; z(z > zmax) = zmax;

scatter(x,y,pointsize,z,'filled')
colormap turbo
colorbar
cleanSubplot
% axis square

end