%% Load data, make LineageTrace object
clear; close all;
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);

% raw data location, modify if not the same as the location of this script
dataDir = scriptPath;
savedir = fullfile(dataDir,'figures');
if ~exist(savedir,'dir'), mkdir(savedir); end

liveDir = 'D:\220928_Smad4GFP_BMP_IWP2_live';
% fixedDir = fullfile(liveDir,'fixed');
fixedDir = fullfile(liveDir,'fixed_aligned\round1');

liveMeta = load(fullfile(liveDir,'meta.mat'));
liveMeta = liveMeta.meta;
livePos = load(fullfile(liveDir,'positions.mat'));
livePos = livePos.positions;

fixedMeta = load(fullfile(fixedDir,'meta.mat'));
fixedMeta = fixedMeta.meta;
fixedPos = load(fullfile(fixedDir,'positions.mat'));
fixedPos = fixedPos.positions;

channelLabels = fixedMeta.channelLabel;

%initialize lineage trace
lt = LineageTrace(livePos, fixedPos, liveMeta, fixedMeta, liveDir, fixedDir);
npos = length(lt.live_position);

%get image size
A = lt.fixed_position(1).loadImage(fixedDir,0,1);
imsize = size(A,[1,2]);
clear A
%number of time points in live imaging
ntime = lt.live_position(1).nTime;
nucChannel = 0;

treatmentTime = 4;
s = strsplit(liveMeta.timeInterval,'min');
tvec = ((0:ntime-1) - treatmentTime)*str2double(s{1})/60;

%% Map live to fixed
close all
positionIdx = 1:npos;
maxDist = 35;
mapPoints(lt, maxDist, struct('alignmentMethod','automatic',...
    'positionIdx',positionIdx));

for ii = positionIdx
    lt.checkAlignment(ii);
    cleanSubplot(18)
    title(sprintf('Position #%d',ii))
    pause
    clf
end
close all

save(fullfile(dataDir,'lt.mat'),'lt')

%% Or load existing LineageTrace object
load(fullfile(dataDir,'lt_correctedcyto.mat'))

%% load tracking results & assign fate
load(fullfile(dataDir,'corrected_tracking','validatedTracking.mat'))
fields = {'XY','NCratio','nucLevel','cytLevel','nucArea','nucZ',...
    'nucMajorAxis','nucMinorAxis','labels'};

hists = cell(1,npos);
for pidx = 1:npos
    ul = UL(pidx,:); ymin = ul(1); xmin = ul(2);
    histories = graphSignalingHistories(lt.live_position(pidx),...
        P(pidx).G,fields,1:length(verified{pidx}));
    start_times = cellfun(@(x) x(1), {histories.Time});
    hists{pidx} = histories(start_times==1);
    
    XYfinal = lt.live_position(pidx).cellData(end).XY;
    mapped = lt.mapped_idxs{pidx};
    nchan = size(lt.fixed_position(pidx).cellData.nucLevel,2);
    
    %only keep histories mapped to fixed data
    cellidxs = cellfun(@(x) x(end), {hists{pidx}.CellIdxs});
    hists{pidx} = hists{pidx}(ismember(cellidxs,mapped));
    
    for ii = 1:length(hists{pidx})
        %mark each cell as verified or not
        hists{pidx}(ii).verified = verified{pidx}(hists{pidx}(ii).CellIdxs(end));
        %determine index of final cell in the history
        xy = hists{pidx}(ii).XY(end,:) + [xmin ymin];
        d = sum((XYfinal - xy).^2,2);
        [~,I] = min(d);
        if I ~= hists{pidx}(ii).CellIdxs(end)
            error('bad indexing')
        end
        fixedIdx = find(mapped == I);
        hists{pidx}(ii).fateMarkers = lt.fixed_position(pidx).cellData.nucLevel(fixedIdx,:);
        hists{pidx}(ii).fixedlabels = lt.fixed_position(pidx).cellData.labels(fixedIdx);
    end
    lt.histories{pidx} = hists{pidx};
end


%% save data matrices as mat file

opts = struct('peakremoval','heuristic','domedfilt',false);
savename = '220928_histories';

%specify the IF channels to save
Mats = lt.histories2mats(1:npos,opts);
vr = cell2mat(cellfun(@(x) [x.verified], lt.histories, 'UniformOutput', false));

save(fullfile(dataDir,[savename,'.mat']),'liveMeta','fixedMeta','Mats','vr')

%% normalize cell fate data
%choose channels for which to save normalized values
channels = [4 3 6 13 14 2 8];
fm = Mats.fateMarkers(channels,:)';
normalizedChannels = fixedMeta.channelLabel(channels);

%load control data
cdata = load(fullfile(liveDir,'fixed_aligned','cntrl','positions'));
cmeta = load(fullfile(liveDir,'fixed_aligned','cntrl','meta'));
cmeta = cmeta.meta;
nc = length(channels);

ncond = length(cmeta.conditions);
ppc = length(cdata.positions)/ncond;
Ms = cell(ncond,1);
for cidx = 1:ncond
    m = cell(ppc,1);
    for cpi = 1:ppc
        pidx = (cidx - 1)*ppc + cpi;
        M = cdata.positions(pidx).cellData.nucLevel(:,channels);
        nanmask = any(isnan(M),2);
        labelmask = cdata.positions(pidx).cellData.labels == 1;
        mask = ~nanmask & labelmask;
        m{cpi} = M(mask,:);
    end
    Ms{cidx} = cell2mat(m);
end

idx1 = 3; idx2 = 4;

fs = 24; lfs = 18;
close all

figure('WindowState','maximized')
np = 100;
pts = NaN(np,nc);
mins = NaN(1,nc); maxes = NaN(1,nc);
normalizedMarkers = NaN(size(fm));

%normalization: would it be better to use averages instead of peaks in the
%density estimate for each channel in each condition?
diffs = NaN(1,nc); mids = NaN(1,nc);
for ii = 1:nc
    m1 = Ms{idx1}(:,ii); m2 = Ms{idx2}(:,ii);
    m1(m1 < 1) = 1; m1 = log(m1 + 1);
    m2(m2 < 1) = 1; m2 = log(m2 + 1);
    mall = [m1;m2];
    xl = axlims([m1;m2],[m1;m2],0.003,false);
    mins(ii) = xl(1); maxes(ii) = xl(2);
    pts(:,ii) = linspace(xl(1),xl(2),np);
    
    subplot_tight(2,4,ii,0.07)
    [f,xi] = ksdensity(m1,pts(:,ii));
    plot(xi,f,'LineWidth',3); hold on
    [~,I] = max(f);
    max1 = xi(I);
    
    [f,xi] = ksdensity(m2,pts(:,ii));
    plot(xi,f,'LineWidth',3); hold off
    [~,I] = max(f);
    max2 = xi(I);
    
    cleanSubplot(fs); axis square
    xlabel(normalizedChannels{ii}); ylabel('density estimate')
%     if ii == nc
%         legend(meta.conditions,'FontSize',lfs,'Location','northeast')
%     end
    
    diffs(ii) = abs(max1 - max2); mids(ii) = mean([max1,max2]);
%     Xnorm(:,ii) = 2*(mall - mids(ii))/diffs(ii);
    v = fm(:,ii); v(v < 1) = 1; v = log(v + 1);
    normalizedMarkers(:,ii) = 2*(v - mids(ii))/diffs(ii);
end

save(fullfile(dataDir,[savename,'.mat']),'liveMeta','fixedMeta',...
    'Mats','vr','channels','normalizedChannels','normalizedMarkers')


%% save as h5 file
%to save:
%   all fields of Mats
%   vr
%   normalizedMarkers
%   channels
%   live channel labels
%   fixed channel labels
%   conditions
%   tvec

% writename = fullfile(dataDir,'220928_histories.h5');
writename = fullfile(dataDir,[savename,'.h5']);

fields = fieldnames(Mats);

%save data from Mats struct
for fi = 1:length(fields)
    field = fields{fi};
    X = Mats.(field);
    h5create(writename, ['/',field], size(X))
    h5write(writename, ['/',field], X)
end

%vector of times in hours
X = tvec;
field = 'tvec';
h5create(writename, ['/',field], size(X))
h5write(writename, ['/',field], X)

%boolean specifying whether each trace has been verified
X = double(vr);
field = 'vr';
h5create(writename, ['/',field], size(X))
h5write(writename, ['/',field], X)

%normalized cell fate markers
X = normalizedMarkers;
field = 'normalizedMarkers';
h5create(writename, ['/',field], size(X))
h5write(writename, ['/',field], X)

%string type fields
%normalized channels labels
X = string(normalizedChannels);
field = 'normalizedChannels';

% h5createwritestr(filename, dataset, str)
h5createwritestr(writename, ['/',field], X);

X = string(fixedMeta.channelLabel);
field = 'fixedChannels';
h5createwritestr(writename, ['/',field], X);

X = string(liveMeta.channelLabel);
field = 'liveChannels';
h5createwritestr(writename, ['/',field], X);

X = string(liveMeta.conditions);
field = 'conditions';
h5createwritestr(writename, ['/',field], X);

%% local function
%needed to write cells of strings to an h5 file
%@author Pavel Komarov pavel@gatech.edu 941-545-7573
function h5createwritestr(filename, dataset, str)

    %"The class of input data must be cellstring instead of char when the
    %HDF5 class is VARIABLE LENGTH H5T_STRING.", but also I don't want to
    %force the user to put braces around single strings, so this.
    if ischar(str)
        str = {str};
    end

    %check whether the specified .h5 exists and either create or open
    %accordingly
    if ~exist(filename, 'file')
        file = H5F.create(filename, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
    else
        file = H5F.open(filename, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
    end

    %set variable length string type
    vlstr_type = H5T.copy('H5T_C_S1');
    H5T.set_size(vlstr_type,'H5T_VARIABLE');

    % There is no way to check whether a dataset exists, so just try to
    % open it, and if that fails, create it.
    try
        dset = H5D.open(file, dataset);
        H5D.set_extent(dset, fliplr(size(str)));
    catch
        %create the intermediate groups one at a time because evidently the
        %API's functions aren't smart enough to be able to do this themselves.
        slashes = strfind(dataset, '/');
        for i = 2:length(slashes)
            url = dataset(1:(slashes(i)-1));%pull out the url of the next level
            try
                H5G.create(file, url, 1024);%1024 "specifies the number of
            catch   %bytes to reserve for the names that will appear in the group"
            end
        end

        %create a dataspace for cellstr
        H5S_UNLIMITED = H5ML.get_constant_value('H5S_UNLIMITED');
%         spacerank = max(1, sum(size(str) > 1));
        spacerank = length(size(str));
        dspace = H5S.create_simple(spacerank, fliplr(size(str)), ones(1, spacerank)*H5S_UNLIMITED);

        %create a dataset plist for chunking. (A dataset can't be unlimited
        %unless the chunk size is defined.)
        plist = H5P.create('H5P_DATASET_CREATE');
        chunksize = ones(1, spacerank);
        chunksize(1) = 2;
        H5P.set_chunk(plist, chunksize);% 2 strings per chunk
        dset = H5D.create(file, dataset, vlstr_type, dspace, plist);

        %close things
        H5P.close(plist);
        H5S.close(dspace);
    end

    %write data
    H5D.write(dset, vlstr_type, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', str);

    %close file & resources
    H5T.close(vlstr_type);
    H5D.close(dset);
    H5F.close(file);
end




