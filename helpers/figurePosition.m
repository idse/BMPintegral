function figpos = figurePosition(width,height)
%make the 4 element position coordinates for a figure of specified width
%and height to display in the center of the screen
%width and height can be specified as separate inputs or as a single
%variable as [width height]

if ~exist('height','var')
    wh = width;
else
    wh = [width,height];
end

ss = get(0,'screensize');
ss = ss(3:4);

figpos = [0.5*(ss - wh),wh];

end