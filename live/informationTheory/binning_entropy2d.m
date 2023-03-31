function H = binning_entropy2d(xy, nbins)
%input should be two-dimesional collection of points of size Nx2 (first
%column is first coordinate, second column is second coordinate)
%nbins gives the number of bins; if nbins is a scalar, this will be the
%number of bins for each coordinate and the total number of 2d bins will be
%nbins^2; alternatively nbins can be specified as a vector with a different
%number of bins for x and y

if length(nbins) == 1
    nx = nbins;
    ny = nbins;
elseif length(nbins) == 2
    nx = nbins(1);
    ny = nbins(2);
else
    error('nbins should have length 1 or 2')
end

xrange = [min(xy(:,1)), max(xy(:,1))];
xdiff = xrange(2) - xrange(1);
xrange = xrange + xdiff*[-0.01, 0.01];
xbins = linspace(xrange(1), xrange(2), nx + 1);

yrange = [min(xy(:,2)), max(xy(:,2))];
ydiff = yrange(2) - yrange(1);
yrange = yrange + ydiff*[-0.01, 0.01];
ybins = linspace(yrange(1), yrange(2), ny + 1);

N = size(xy,1);

[P,~,~] = histcounts2(xy(:,1), xy(:,2), xbins, ybins);

P = P(:);
P = P(P > 0)/N;

H = -sum(P.*log2(P));




end