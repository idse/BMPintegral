function [xl, yl] = axlims(xdata,ydata,tol,setnow)

%default tolerance (set to zero to make axes tight around data without
%cropping anything out)
if ~exist('tol','var')
    tol = 0;
end

if ~exist('setnow','var')
    setnow = true;
end

%can specify different tolerances in x and y if desired
if numel(tol) == 2
    tolx = tol(1); toly = tol(2);
elseif numel(tol) == 1
    tolx = tol; toly = tol;
end

n = length(xdata);
%probably not the most computationally efficient way to do this
xs = sort(xdata);
ys = sort(ydata);

nx = n - sum(isnan(xs)); %don't include NaN values in calculation
lhx = [max(1,floor(nx*tolx)), ceil(nx*(1-tolx))];
xl = [xs(lhx(1)),xs(lhx(2))];

ny = n - sum(isnan(ys)); %don't include NaN values in calculation
lhy = [max(1,floor(ny*toly)), ceil(ny*(1-toly))];
yl = [ys(lhy(1)),ys(lhy(2))];

if setnow
    set(gca,'xlim',xl)
    set(gca,'ylim',yl)
end

end
