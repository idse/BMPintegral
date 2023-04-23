function [X,Err,Rs,edges] = newRadialProfile(positions,chan,opts)

nc = length(chan);

if isfield(opts,'fields')
    fields = opts.fields;
else
    fields = repmat({'nucLevel'},1,nc);
end

if isfield(opts,'t')
    t = opts.t;
else
    t = positions(1).nTime;
end

if isfield(opts,'normalizeToDAPI')
    DAPInorm = opts.normalizeToDAPI;
    if numel(DAPInorm) == 1 && nc > 1
        DAPInorm = repmat(DAPInorm,1,nc);
    end
else
    DAPInorm = false(1,nc);
end

if isfield(opts,'nbin')
    nbin = opts.nbin;
else
    nbin = 35;
end

if ~isfield(opts,'normalize')
    opts.normalize = false;
end

if ~isfield(opts,'suppressoutput')
    opts.suppressoutput = false;
end

if ~isfield(opts,'maxval')
    maxval = Inf(1,nc);
else
    maxval = opts.maxval;
end
if ~isfield(opts,'minval')
    minval = -Inf(1,nc);
else
    minval = opts.minval;
end

if ~isfield(opts,'avgmethod')
    opts.avgmethod = 'median';
end
%percentage overlap in radial bins (overlapping bins makes the profile vary
%smoothly)
% if ~isfield(opts,'overlap') || isempty(opts.overlap)
%     overlap = 0;
% else
%     overlap = opts.overlap;
% end

%make the number of angular bins an input?
%for now, use 50 as the default
theta = linspace(0,2*pi,50);
nth = length(theta);

nucChannel = positions(1).nucChannel + 1;
rad = positions(1).radiusMicron;
xres = positions(1).radiusMicron/positions(1).radiusPixel;

npos = length(positions);
data = cell(npos,1); D = cell(npos,1);

for pidx = 1:npos    
    %get centered xy positions of all cell centroids in microns
    XY = positions(pidx).cellData(t).XY;
    cm = positions(pidx).center(t,:);
    xy = (XY - cm)*xres;
    
    %convert to polar coordinates (use different inverse trig functions in
    %different quadrants)
    R = sqrt(sum(xy.^2,2));
    Th = NaN(size(R));
    idx = find(xy(:,2) >= 0);
    Th(idx) = acos(xy(idx,1)./R(idx));
    idx = find(xy(:,1) <= 0 & xy(:,2) <= 0);
    Th(idx) = pi + atan(xy(idx,2)./xy(idx,1));
    idx = find(xy(:,1) >= 0 & xy(:,2) <= 0);
    Th(idx) = 2*pi + atan(xy(idx,2)./xy(idx,1));
    
    %find max radial distance less than nominal among cells within angular
    %bins
    xedge = NaN(nth-1,1); yedge = xedge;
    for ii = 1:nth-1
        idx = Th > theta(ii) & Th < theta(ii+1) & R < rad;
        rtemp = R(idx); ttemp = Th(idx);
        [~,I] = max(rtemp);
        xedge(ii) = rtemp(I)*cos(ttemp(I));
        yedge(ii) = rtemp(I)*sin(ttemp(I));
    end
    %interpolate between edge cells for more accurate edge distance
    %calculation (closer to a continuous boundary around the colony)
    Xedge = cell(nth-1,1); Yedge = cell(nth-1,1);
    for ii = 1:nth-1
        x1 = xedge(ii); y1 = yedge(ii);
        if ii < nth - 1
            x2 = xedge(ii+1); y2 = yedge(ii+1);
        elseif ii == nth - 1
            x2 = xedge(1); y2 = yedge(1);
        end
        xs = linspace(x1,x2,100); ys = linspace(y1,y2,100);
        Xedge{ii} = xs(:); Yedge{ii} = ys(:);
    end
    xedge = cell2mat(Xedge); yedge = cell2mat(Yedge);
    
    %find distance to nearest point on calculated colony boundary
    [~, dist] = knnsearch([xedge,yedge],xy);
    D{pidx} = dist(:);
    
    %collect data for each channel
    data{pidx} = NaN(size(XY,1),nc);
    DAPI = positions(pidx).cellData(t).nucLevel(:,nucChannel);
    for ci = 1:nc
        vals = positions(pidx).cellData(t).(fields{ci})(:,chan(ci));
        if DAPInorm(ci)
            data{pidx}(:,ci) = vals./DAPI;
        else
            data{pidx}(:,ci) = vals;
        end
    end
end

dist = cell2mat(D);
Rmax = max(dist);
data = cell2mat(data);
ncells = length(dist); %number of cells
cpb = ncells/nbin; %cells per bin
% opb = round(overlap*cpb); %number of overlapping cells per bin (on each side)

if ~opts.suppressoutput
    fprintf('ncells = %d, cells per bin = %d\n',ncells,round(cpb))
end

[sorted,I] = sort(dist);
data = data(I,:);

Rs = NaN(nbin,1); X = NaN(nbin,nc); Err = NaN(nbin,nc);
edges = NaN(nbin+1,1);
for bi = 1:nbin
%     ind1 = max(1,round((bi-1)*cpb) + 1 - opb);
%     ind2 = min(ncells,round((bi-1)*cpb)+1 + opb);
%     inds = ind1:ind2;
%     inds = round((bi-1)*cpb) + 1:round((bi-1)*cpb)+1;
    inds = round((bi-1)*cpb) + 1:round(bi*cpb);
    Rs(bi) = mean(sorted(inds));
    edges(bi) = sorted(inds(1));
    for ci = 1:nc
        vals = data(inds,ci);
        vals(vals > maxval(ci) | vals < minval(ci)) = NaN;
        
        if strcmp(opts.avgmethod,'mean')
            X(bi,ci) = mean(vals,'omitnan');
        elseif strcmp(opts.avgmethod,'median')
            X(bi,ci) = median(vals,'omitnan');
        else
            error('avgmethod should be mean or median')
        end
        
        Err(bi,ci) = std(vals,'omitnan');
    end
end
edges(end) = Rmax;

if opts.normalize
    for ci = 1:nc
        limits = [min(X(:,ci)),max(X(:,ci))];
        X(:,ci) = (X(:,ci) - limits(1))./(limits(2) - limits(1));
        Err(:,ci) = Err(:,ci)./(limits(2) - limits(1));
    end
end







end