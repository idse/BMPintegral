function [imgs, clabel, savename] = profile2pattern(V,rs,opts)

if ~iscell(V)
    V = {V};
end

ncond = length(V);
nc = size(V{1},2);

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
if isfield(opts,'colorMode')
    colorMode = opts.colorMode;
else
    colorMode = 'rgb';
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

if ~isfield(opts,'order')
    order = zeros(1,3);
    order(1:min(nc,3)) = 1:min(nc,3);
else
    order = opts.order;
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
        ci = order(jj);
        im = B;
        if ci > 0
            vals = V{ii}(I,ci); %get expression values in channel ci in the same order as r
            for bi = 1:length(edges)-1
                mask = (DIST < edges(bi)) & (DIST >= edges(bi+1));
                im(mask) = vals(bi);
            end
        end
        img(:,:,jj) = im;
    end
%     img = imgs{ii};
    if strcmpi(colorMode,'cmy')
        K = min(1 - img,[],3); %black level
        img = (1 - img).*(1 - K); %convert rgb to cmy
    end
    for jj = 1:3
        im = img(:,:,jj);
        im(bmask) = 1;
        img(:,:,jj) = im;
    end
    imgs{ii} = img;
end

% error('troubleshooting')

if strcmpi(colorMode,'rgb')
    cols = {'{red}','{green}','{blue}'};
    savecols = {'r','g','b'};
elseif strcmpi(colorMode,'cmy')
    cols = {'{cyan}','{magenta}','{yellow}'};
    savecols = {'c','m','y'};
end
chans = channels(order(order > 0));
cols = cols(order > 0);
savecols = cell2mat(savecols(order > 0));

%handle case in which we combine one rgb color with one cmy color; assume
%that there is never more than one empty color
if order(1) == order(2)
    cols = {'{yellow}','{blue}'};
    savecols = 'yb';
    chans = channels(order([1,3]));
elseif order(2) == order(3)
    cols = {'{red}','{cyan}'};
    savecols = 'rc';
    chans = channels(order([1,2]));
elseif order(3) == order(1)
    cols = {'{green}','{magenta}'};
    savecols = 'gm';
    chans = channels(order([2,3]));
end

clabel = cell2mat(strcat('\color',cols,chans," "));
clabel = clabel(1:end-1);

% savecol = strrep(strrep(cell2mat(cols),'}',''),'{','');
savename = strcat('_',savecols,cell2mat(strcat('_',chans)));










end