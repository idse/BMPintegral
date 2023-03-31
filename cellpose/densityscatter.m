function densityscatter(x,y,tol)

if ~exist('tol','var')
    tol = 0.02;
end

[f,xi] = ksdensity([x,y]);
np = sqrt(size(xi,1));
Xi = reshape(xi(:,1),[np,np]); Yi = reshape(xi(:,2),[np,np]); F = reshape(f,[np,np]);
z = interp2(Xi,Yi,F,x,y);


n = length(z);
zs = sort(z);
zmin = zs(max(ceil(n*tol),1));
zmax = zs(floor(n*(1-tol)));
z(z < zmin) = zmin; z(z > zmax) = zmax;

scatter(x,y,15,z,'filled')
cleanSubplot
axis square

end