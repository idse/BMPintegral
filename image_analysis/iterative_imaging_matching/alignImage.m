function aligned = alignImage(im1, im2, shiftyx)
%after mapping points, overlay the (shifted) images and show
%the link between fixed and live points

m1 = size(im1,1); n1 = size(im1,2); m2 = size(im2,1); n2 = size(im2,2);
mn = min([m1,n1], [m2,n2]);
% m = mn(1); n = mn(2);
% im2 = im2(1:mn(1), 1:mn(2));

dy = shiftyx(1);
dx = shiftyx(2);

% xinrange = max(1,1+dx):min(m,m+dx);
% xoutrange = max(1,1-dx):min(m,m-dx);
% yinrange = max(1,1+dy):min(n,n+dy);
% youtrange = max(1,1-dy):min(n,n-dy);

xin1 = max(1,1+dx);
xin2 = min(m2,m1+dx);
xdiff = xin2 - xin1 + 1;
xinrange = xin1:xin2;

yin1 = max(1,1+dy); yin2 = min(n2,n1+dy);
ydiff = yin2 - yin1 + 1;
yinrange = yin1:yin2;

xout1 = max(1,1-dx); xout2 = min(m1,xout1 + xdiff - 1);
xoutrange = xout1:xout2;

yout1 = max(1,1-dy); yout2 = min(n1,yout1 + ydiff - 1);
youtrange = yout1:yout2;

aligned = zeros(size(im1),class(im2));
aligned(xoutrange, youtrange) = im2(xinrange, yinrange);

end