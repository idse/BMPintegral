function v = cleanHistory_v2(histories,idx,opts)
%function to build a single-cell history of a specific cell in a given
%channel for a desired field
%first, NaN values are interpolated
%then a heuristic function is used to identify cell divisions and the peaks
%due to cell division are removed and values are interpolated
%finally, the history is smoothed with a first-order median filter
%%% inputs %%%
%histories is a struct of single-cell histories
%idx is the index of the desired cell
%field is the cellData field of interest
%channel is the cellData channel (applicable if field is NCratio, nucLevel,
%or cytLevel)

if ~exist('opts','var')
    opts = struct;
end

%if a number is specified instead of a struct for opts, that number is the
%channel, all other options are set to default
if ~isa(opts,'struct')
    channel = opts;
    opts = struct;
else
    channel = opts.channel;
end


if ~isfield(opts,'field')
    opts.field = 'NCratio';
end
if ~isfield(opts,'nucChannel')
    nucChannel = 0;
else
    nucChannel = opts.nucChannel;
end
if ~isfield(opts,'domedfilt')
    opts.domedfilt = false;
end
if ~isfield(opts,'peakremoval')
    opts.peakremoval = 'heuristic';
end


minval = 0;

if strcmp(opts.field,'NCratio') && ~(channel==nucChannel)
    maxval = 2;
else
    maxval = Inf;
end

time = histories(idx).Time;
ntime = time(end);
tvec = (1:ntime)';

%collect nuclear statistics and interpolate missing data
A = interp1(time,histories(idx).nucArea,tvec);
ratio = histories(idx).nucMajorAxis./histories(idx).nucMinorAxis;
R = interp1(time,ratio,tvec);
N = interp1(time,histories(idx).nucLevel(:,nucChannel+1),tvec);
X = interp1(time,histories(idx).(opts.field)(:,channel+1),tvec);

if strcmp(opts.peakremoval,'heuristic')
    %heuristic cell division detection
    criteria = (A<mean(A,'omitnan') - 2*std(A,'omitnan')) + (R>2.5) +...
        (X>mean(X,'omitnan') + 2*std(X,'omitnan')) +...
        (N>mean(N,'omitnan') + 2*std(N,'omitnan'));
    criteria = (criteria > 1) | X > maxval | X < minval;
    criteria = conv(criteria,ones(8,1),'same');
    %remove values around cell division peaks and interpolate
    X = X(criteria == 0);
    vt = tvec(criteria == 0);
    v = interp1(vt,X,tvec,'linear');
elseif strcmp(opts.peakremoval,'labelbased')
    %identify cell division based on labels in celldata (generally based on
    %ilastik object classification of cells using image data)
    L = ones(ntime,1);
    L(time) = histories(idx).labels;
    criteria = L > 1;
    
    criteria = conv(criteria,ones(8,1),'same');
    X = X(criteria == 0);
    vt = tvec(criteria == 0);
    v = interp1(vt,X,tvec,'linear');
elseif strcmp(opts.peakremoval,'none')
    %do not remove cell division peaks
    v = X;
end
%interpolate NaN values
nonans = find(~isnan(v));
v(nonans(end):end) = v(nonans(end)); %extrapolate (nearest, not linear)
v(1:nonans(1)) = v(nonans(1)); %extrapolate (nearest, not linear)

if opts.domedfilt
    v = medfilt1(v,'omitnan','truncate');
end

%additional NaN value interpolation -> why do we need this? is it possible
%to get nans from median filtering if there were none before, or is this
%not needed anymore?
if any(isnan(v))
    nonans = find(~isnan(v));
    v1 = v(nonans);
    v = interp1(nonans,v1,tvec,'linear');
    v(1:nonans(1)) = v(nonans(1));
    v(nonans(end):end) = v(nonans(end));
end



end

