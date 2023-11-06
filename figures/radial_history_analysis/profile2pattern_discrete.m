function [imgs, savename] = profile2pattern_discrete(V,rs,opts)

if ~iscell(V)
    V = {V};
end

ncond = length(V);
uv = unique(cell2mat(V));
nc = length(uv);

if ~exist('opts','var')
    opts = struct;
end

if ~isfield(opts,'Rmax')
    Rmax = 350;
else
    Rmax = opts.Rmax;
end

%allow specification of either rgb or cmy colors; should also be able to
%handle red/cyan or green/magenta two-color combinations, by using rgb mode
%and letting one of the channels occupy both blue and green or both blue
%and red colors
if isfield(opts,'colors')
    colors = lines(nc);
else
    colors = opts.colors;
end

if ~isfield(opts,'channels')
    channels = cellstr(num2str((1:nc)'))';
else
    channels = opts.channels;
end

if ~iscell(rs)
    r = rs;
    rs = cell(size(V));
    for ii = 1:ncond
        rs{ii} = r;
    end
end

imgs = cell(size(V));
imsize = round(2.3*Rmax)*ones(1,2);
B = zeros(imsize);
DIST = zeros(imsize);
DIST(round(imsize(1)/2),round(imsize(2)/2)) = 1;
DIST = bwdist(DIST);
%background mask -> make the image inside a black circle surrounded by white
bmask = DIST >= 1.05*Rmax;%imsize(1)/2;

for ii = 1:ncond
    r = rs{ii}(:);
    [r,I] = sort(r); %put r in ascending order
    if strcmp(opts.rmode,'edges')
        rmax = max(r);
        edges = Rmax/rmax*(max(r) - r);
        I = I(1:end-1);
    elseif strcmp(opts.rmode,'radius')
        edges = Rmax - [0; 0.5*(r(2:end) + r(1:end-1)); Rmax];
    end
%     edges = Rmax - [0; 0.5*(r(2:end) + r(1:end-1)); Rmax];
    img = zeros(imsize(1),imsize(2),3);
    for jj = 1:3
        im = B;
        vals = V{ii}(I); %get expression values in channel ci in the same order as r
        for bi = 1:length(edges)-1
            cval = colors(uv == vals(bi),jj); %color value
            mask = (DIST < edges(bi)) & (DIST >= edges(bi+1));
            im(mask) = cval;
        end
        img(:,:,jj) = im;
    end
%     img = imgs{ii};
    for jj = 1:3
        im = img(:,:,jj);
        im(bmask) = 1;
        img(:,:,jj) = im;
    end
    imgs{ii} = img;
end

% error('troubleshooting')

% clabel = cell2mat(strcat('\color',cols,chans," "));
% clabel = clabel(1:end-1);

% savecol = strrep(strrep(cell2mat(cols),'}',''),'{','');
savename = strcat('_',cell2mat(strcat('_',channels)));










end