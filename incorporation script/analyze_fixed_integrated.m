clear; close all;

scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);

% raw data location, modify if not the same as the location of this script
dataDir = scriptPath;
cd(dataDir);

% determine MP or Disordered
imageType = 'MP'; % 'MP' or 'Disordered'

manMeta = struct();
manMeta.nucChannel = 0;
manMeta.channelLabel = {'DAPI','BRDU','BRA','YAP'};

manMeta.conditions = {'80K','160K','240K'};

if strcmp(imageType,'MP')
    radii = [700 700 700]/2;
end

manMeta.posPerCondition = 4;
manMeta.nWells = numel(manMeta.conditions);
manMeta.nPositions = manMeta.posPerCondition*manMeta.nWells;

JOB_Nikon = false;
multiposition_Dragonfly = false;

if JOB_Nikon == 1 || multiposition_Dragonfly == 1
    filelistAll = {};
    if JOB_Nikon == 1
        filelistAll{1} = dir(fullfile(dataDir,'*.nd2')); % read all NIKON JOB files
    elseif multiposition_Dragonfly == 1
        filelistAll{1} = dir(fullfile(dataDir,'*.ims')); % read all Dragonfly multiposition files
    end

    filelist = {};
    for i = 1:numel(manMeta.conditions)
        filelist{i} = filelistAll{1}((i-1)*manMeta.posPerCondition+1:i*manMeta.posPerCondition);
    end
else % determine filelist condition by condition
    filelist = {};
    filelist{1} = dir(fullfile(dataDir,'*80K*.ims'));
    filelist{2} = dir(fullfile(dataDir,'*160K*.ims'));
    filelist{3} = dir(fullfile(dataDir,'*240K*.ims'));
end

coloniesPerSize = cellfun(@numel,filelist);
posbins = [0 cumsum(coloniesPerSize)];

posfile = fullfile(dataDir,['positions' [manMeta.channelLabel{:} '.mat']]);
if exist(posfile)
    load(posfile, 'positions');
    disp(['loading ' posfile]);
end

statsfile = fullfile(dataDir,['stats' [manMeta.channelLabel{:} '.mat']]);
if exist(statsfile)
    load(statsfile, 'stats');
    disp(['loading ' statsfile]);
end

countfile = fullfile(dataDir,['counts' [manMeta.channelLabel{:} '.mat']]);
if exist(countfile)
    load(countfile, 'counts');
    disp(['loading ' countfile]);
end

metafile = fullfile(dataDir,['meta'  [manMeta.channelLabel{:} '.mat']]);
if exist(metafile)
    load(metafile,'meta');
    disp(['loading ' metafile]);
else
    fname = filelist{1}(1).name;
    subdir = filelist{1}(1).folder;
    meta = Metadata(fullfile(subdir,fname), manMeta);
    save(metafile, 'meta');
end

%% make xsections of different channels
pList = [1,5,9];

count = 0; % use first position to define lim
for pi = pList
    count = count + 1;
    fname = filelist{floor(pi/manMeta.posPerCondition)+1}(mod(pi,manMeta.posPerCondition)).name;
    img = readStack(fullfile(dataDir,fname));

    yinds = {1250:1260; 1000:1010}; %change yinds if preferred
    for i = 1:numel(yinds)
        RGBset = [2 3 4];
        nucChannel = 1;
        tol = [0.01 0.99; 0.01 0.99; 0.01 0.99; 0.1 0.99];
        zcorrection = true;
        if count == 1
            disp('setting Ilim');
            Ilim = makeXsection(img, fname, meta, yinds{i}, RGBset, nucChannel, tol,zcorrection);
        else
            makeXsection(img, fname, meta, yinds{i}, RGBset, nucChannel, tol,zcorrection, Ilim);
        end
    end
end

%% extractData
opts = struct('cytoSize',4,'cytoMargin',2,'IoU',0.5,'Zmethod','new');

if strcmp(imageType,'MP')
    positions(manMeta.nPositions) = Colony();
elseif strcmp(imageType,'Disordered')
    positions(manMeta.nPositions) = Position();
end

% for all conditions (aka wells)
tic
for ci = 1:manMeta.nWells

    % for all positions belonging to that condition
    for cpi = 1:manMeta.posPerCondition

        pi = manMeta.posPerCondition*(ci-1) + cpi;

        fname = filelist{ci}(cpi).name;
        subdir = filelist{ci}(cpi).folder;

        meta = Metadata(fullfile(subdir,fname), manMeta);

        meta.nZslices = meta.nZslices;
        % -1 is only because fusion stitcher attaches a black z-slice at the end
        if contains(fname,'FusionStitcher')
            meta.nZslices = meta.nZslices - 1;
        end

        %if there is a stitched 4D tiff with a name matching this one, use
        %that instead of the FusionStitcher file
        tifname = strrep(strrep(fname,'FusionStitcher','Stitched'),'.ims','.tif');
        if exist(fullfile(subdir,tifname),'file') == 2
            fname = tifname;
            r = bfGetReader(fullfile(subdir,tifname));
            %update metadata
            meta.nZslices = r.getSizeZ(); %don't subtract 1 for stitched tiff stacks
            meta.xSize = r.getSizeX();
            meta.ySize = r.getSizeY();
            fnameformat = meta.filenameFormat;
            meta.filenameFormat =...
                strrep(strrep(fnameformat,'FusionStitcher','Stitched'),'.ims','.tif');
        end

        % id ~= pi in this case
        if ~isempty(regexp(fname, '_F[0-9]','Once'))
            extension = fname(end-3:end);
            s = strsplit(fname,{'_F',extension});
            id = uint16(str2double(s{end-1})) + 1;
        else
            id = pi;
        end
        fprintf('Position %d\n',pi)
        disp(fname)
        if strcmp(imageType,'MP')
            positions(pi) = Colony(meta, id);
            positions(pi).setRadius(radii(ci), meta.xres);
            positions(pi).well = ci;
            positions(pi).dataChannels = 1:positions(pi).nChannels;

            positions(pi).cellData = extractDataMultiZ(subdir, fname, meta, opts);
            positions(pi).setCenter();
            positions(pi).makeRadialAvgSeg();
            positions(pi).ncells = size(positions(pi).cellData.XY,1);
        elseif strcmp(imageType,'Disordered')
            positions(pi) = Position(meta, id);
            positions(pi).dataChannels = 1:positions(pi).nChannels;

            positions(pi).cellData = extractDataMultiZ(subdir, fname, meta, opts);
            positions(pi).ncells = size(positions(pi).cellData.XY,1);
        end
    end
end
toc
save(posfile,'positions');

%% make xsections of segmentation
%save 1 position per condition
tic
for ci = 1:manMeta.nWells
    cpi = 1; % No.cpi colony in each condition
    pi = manMeta.posPerCondition*(ci-1) + cpi;
    fname = positions(pi).filename;
    prefix = fname(1:end-4);
    extension = fname(end-3:end);
    img = readStack(fullfile(dataDir,fname));

    LMfname = fullfile(dataDir,[prefix '_LabelMatrix3D.mat']);
    LM = load(LMfname);
    LM = LM.LM;

    s = strsplit(fname, {'_Stitched','_FusionStitcher',extension});
    meta = load(fullfile(dataDir,[[s{:}] '_zslices'],'meta.mat'));
    meta = meta.meta;
    nuclearim = squeeze(img(:,:,meta.nucChannel+1,1:meta.nZslices));
    %change yinds and xinds if preferred
    yinds = 1400;
    xinds = 1:size(img,1);
    makeSegmentationCrossSection(nuclearim,LM,meta,yinds,xinds)
    saveas(gcf, fullfile(dataDir, [prefix  'segmentationXsection_y' num2str(yinds) '.png']));
    figure()
    MIP = max(nuclearim,[],3);
    imshow(imadjust(MIP,stitchedlim(MIP)));
    x = [min(yinds),max(yinds)];
    y = [min(xinds),max(xinds)];
    line(x,y,'Color','r','LineWidth',2)
    saveas(gcf, fullfile(dataDir, [prefix  'segmentationXsection_y_overview' num2str(yinds) '.png']));
end
toc

%% check if you extract position Intensity right

piSerious = [1,5,9];
ci = 1; % define nucChannel you want to check
checkOPT = 'nuc'; % choose from 'nuc' 'cyto' 'NCratio'
for pi = piSerious
    condi = ceil(pi / meta.posPerCondition);

    tol = 0.01;
    img = max(positions(pi).loadImage(dataDir,ci-1,1),[],3);

    % saturate ALL PIXEL according to tol
    imgGREY = imadjust(img,stretchlim(img,[tol,1-tol]),[]);
    imgBLACK = imbinarize(img, 1);

    nucLevel = positions(pi).cellData.nucLevel;
    cytLevel = positions(pi).cellData.cytLevel;
    background = positions(pi).cellData.background;
    XY = positions(pi).cellData.XY;

    cNuc = nucLevel(:,ci)-background(ci);
    cCyt = cytLevel(:,ci)-background(ci);
    NC = cNuc./cCyt;

    if strcmp(checkOPT,'nuc')
        z = cNuc;
    elseif strcmp(checkOPT,'cyto')
        z = cCyt;
    elseif strcmp(checkOPT,'NCratio')
        z = NC;
    end

    % saturate ALL CELLS according to tol
    n = length(z);
    zs = sort(z);
    zmin = zs(max(ceil(n*tol),1));
    zmax = zs(floor(n*(1-tol)));
    z(z < zmin) = zmin; z(z > zmax) = zmax;
    z = (z-min(z))/(max(z)-min(z));

    figure()
    imshow(imgGREY)
    hold on
    s = scatter(XY(:,1), XY(:,2), 20, z*[1, 1, 1], "filled");
    hold off
    %s.MarkerEdgeColor = 'm';
    title(['Intensity in ',meta.channelLabel{ci}])
    saveas(gcf,['Intensity_check_' meta.channelLabel{ci} '_' manMeta.conditions{condi} '_Pos' num2str(pi) '.jpg']);
end

%% CLEAN UP positions & VISUALIZE cells
checkPos = true; % set if you want to VISUALIZE cells
posCLEANUP = false; % set if you want to do clean up

badidxPos = cell(numel(positions),1);
badidxCon = cell(numel(manMeta.conditions),1);
for pos= 1:numel(positions)
    condi = ceil(pos/meta.posPerCondition);
    i = 1; % define which channel you do CLEAN UP according to

    background = positions(pi).cellData.background;
    nucLevel = positions(pos).cellData.nucLevel(:,i)-background(i);
    nucArea = positions(pos).cellData.nucArea(:,1);

    stdFOLD = 2;
    % throwout brightest DAPI (dying/dead cells)
    badidxPos{pos} = (nucLevel > mean(nucLevel) + stdFOLD*std(nucLevel));
    % throw out darkest DAPI (usually segmentation error)
    badidxPos{pos} = badidxPos{pos} | (nucLevel < mean(nucLevel) - stdFOLD*std(nucLevel));
    % throw out cells with too small nucArea, also can use nucAxis, nucCircularity etc
    badidxPos{pos} = badidxPos{pos} | (nucArea < mean(nucArea) - stdFOLD*std(nucArea));
    if strcmp(imageType,'MP')
        % throw out cells outside the colony edges
        badidxPos{pos} = badidxPos{pos}; (sqrt(sum((positions(pos).cellData.XY-positions(pos).center).^2,2))*meta.xres > positions(pos).radiusMicron);
    end
    if posCLEANUP == 0
        badidxCon{condi} = [badidxCon{condi};badidxPos{pos}];
    end

    if checkPos == 1 % check how CLEAN UP works
        figure()
        img = max(positions(pos).loadImage(dataDir,i-1,1),[],3);
        img = mat2gray(img, [150 4000]);
        imshow(img)
        hold on
        XY = positions(pos).cellData.XY;
        XYbad = positions(pos).cellData.XY(badidxPos{pos},:);
        scatter(XY(:,1),XY(:,2), 7,'filled','c')
        scatter(XYbad(:,1),XYbad(:,2), 7,'filled','m')
        hold off
        saveas(gcf, fullfile(dataDir, ['CLEAN UP position' num2str(pos) '.png']));
    end

    if posCLEANUP == 1 % throw out badIdx cells
        positions(pos).ncells = positions(pos).ncells - sum(badidxPos{pos});
        positions(pos).cellData.cytMaskArea = positions(pos).cellData.cytMaskArea(~badidxPos{pos},:);
        positions(pos).cellData.cytLevel = positions(pos).cellData.cytLevel(~badidxPos{pos},:);
        positions(pos).cellData.NCratio = positions(pos).cellData.NCratio(~badidxPos{pos},:);
        positions(pos).cellData.XY = positions(pos).cellData.XY(~badidxPos{pos},:);
        positions(pos).cellData.nucArea = positions(pos).cellData.nucArea(~badidxPos{pos},:);
        positions(pos).cellData.nucOrientation = positions(pos).cellData.nucOrientation(~badidxPos{pos},:);
        positions(pos).cellData.nucMajorAxis = positions(pos).cellData.nucMajorAxis(~badidxPos{pos},:);
        positions(pos).cellData.nucMinorAxis = positions(pos).cellData.nucMinorAxis(~badidxPos{pos},:);
        positions(pos).cellData.nucCircularity = positions(pos).cellData.nucCircularity(~badidxPos{pos},:);
        positions(pos).cellData.nucZ = positions(pos).cellData.nucZ(~badidxPos{pos},:);
        positions(pos).cellData.nucLevel = positions(pos).cellData.nucLevel(~badidxPos{pos},:);
        badidxCon{condi} = [badidxCon{condi};badidxPos{pos}(~badidxPos{pos},:)];
    end

end
%% get stats
stats = cellStats(positions, meta, positions(1).dataChannels);
% automaticly get thresholds
confidence = 0.95;
conditions = 1:numel(meta.conditions);
whichthreshold = []; %[1 1 2]
stats.getThresholds(confidence, conditions, whichthreshold);

% Make nucHistogram channel by channel
conditionIdx = 1:numel(meta.conditions); % [1,3,5] Specify conditions for Histograms
tolerance = 0.01;
nbins = 50;
stats.makeHistograms(nbins, tolerance);
for channelIdx = 1:numel(manMeta.channelLabel)
    options = struct('channelIndex',channelIdx, 'cytoplasmic', false,...
        'cumulative',false, 'time', 1,...
        'conditionIdx',conditionIdx,...
        'titlestr',meta.channelLabel{channelIdx}...
        );
    figure()
    stats.plotDistributionComparison(options)
    saveas(gcf, fullfile(dataDir, [meta.channelLabel{channelIdx(1)} '_dist.png']));
end
%ylim([0 0.02])
%xlim([0 1000])

save(statsfile, 'stats');
%stats.exportCSV(meta); export to CSV
%% Create combined nucHistogram channel by channel, If mannually thresholding, use this as a referrence
nucLevelCombined = [];
XYCombined = [];
sampleIdx = [];
for i = 1:numel(manMeta.conditions)
    nucLevel = stats.nucLevel{i};
    nucLevelCombined = [nucLevelCombined; nucLevel];
end
for i=1:numel(manMeta.channelLabel)
    figure()
    histogram(nucLevelCombined(:,i));
    xlabel('Intensity');
    ylabel('Population');
    title([manMeta.channelLabel{i},' All'])
    %     xlim([0,3000]);
    %     ylim([0,3000]);
    saveas(gcf, fullfile(dataDir, [meta.channelLabel{i} '_tot_dist.jpg']));
end

%% SAMPLE QC PLOT
% check that scatterplots of different samples in same conditions look
% similar

condi = 1; % define the condition you want to check
i = 2; % define two channels you want to compare
j = 3; % define two channels you want to compare

X = stats.nucLevel{condi}(:,i);
Y = stats.nucLevel{condi}(:,j);
X = log(1+X/mean(X));
Y = log(1+Y/mean(Y));
scatter(X, Y, 10, stats.sample{condi}, 'filled')
% xlim([0 1.5]);
% ylim([0 3]);

%% make intensity radial profile (MP Only)
% 1. set 'nucChannels' 'cytChannels' 'ratioChannels' to plot different value
% 2. default interpmethod is 'linear', method 'nearest' is for sparse colony center
% 3. stdMethod: 'perColony' 'neighborCells', usually perColony give less STD
options = struct('ratioChannels', 1:3, 'ratioMode', 'C:N',...
    'normalize', false, 'std', true, 'FontSize', 15, 'legend',true,...
    'pointsPerBin', 200, 'interpmethod', 'nearest', 'stdMethod', 'perColony');
res = {};
% first round without normalization (to collect overall max and min of all), keep normalize = false
for condi = 1:numel(meta.conditions)
    clf
    options.conditionIdx = condi;
    options.colonyRadius = radii(condi);
    res{condi} = plotRadialProfiles(stats, meta, options);
    if condi == 1
        maxvalsNucAll = max(res{condi}.nuc_profile);
        minvalsNucAll = min(res{condi}.nuc_profile);

        maxvalsCytAll = max(res{condi}.cyt_profile);
        minvalsCytAll = min(res{condi}.cyt_profile);

        maxvalsRatioAll = max(res{condi}.ratio_profile);
        minvalsRatioAll = min(res{condi}.ratio_profile);
    else
        maxvalsNucAll = max(max(res{condi}.nuc_profile), maxvalsNucAll);
        minvalsNucAll = min(min(res{condi}.nuc_profile), minvalsNucAll);

        maxvalsCytAll = max(max(res{condi}.cyt_profile), maxvalsCytAll);
        minvalsCytAll = min(min(res{condi}.cyt_profile), minvalsCytAll);

        maxvalsRatioAll = max(max(res{condi}.ratio_profile), maxvalsRatioAll);
        minvalsRatioAll = min(min(res{condi}.ratio_profile), minvalsRatioAll);
    end
end
close all
% second round with normalization
nuclimits = [minvalsNucAll' maxvalsNucAll'];
options.nuclimits = nuclimits;

cytlimits = [minvalsCytAll' maxvalsCytAll'];
options.cytlimits = cytlimits;

ratiolimits = [minvalsRatioAll' maxvalsRatioAll'];
options.ratiolimits = ratiolimits;

options.normalize = true; % set normalization here

for condi = 1:numel(meta.conditions)
    options.conditionIdx = condi;
    options.colonyRadius = radii(condi);
    figure()
    res{condi} = plotRadialProfiles(stats, meta, options);
    ylim([0 1.1]);
    xlim([0 radii(condi)]);
    axis square
    title(meta.conditions{condi})
    saveas(gcf, ['radialProfile_' meta.channelLabel{:} '_' meta.conditions{condi} '.png'])
end
%% make intensity radial profile, one colony (MP Only)
% 1. set 'nucChannels' 'cytChannels' 'ratioChannels' to plot different value
% 2. default interpmethod is 'linear', 'nearest' is for sparse colony center
% 3. stdMethod: 'perColony' 'neighborCells', usually perColony give less STD
% 4. decrease pointsPerBin when you analysis single colony
posList = [1,5,9];
norMethodPos = 1; % 1 - normalize individually; 2 - normalize according to all position
options = struct('ratioChannels', 1:3, 'ratioMode', 'C:N',...
    'normalize', false, 'std', true, 'FontSize', 15, 'legend',true,...
    'pointsPerBin', 100, 'interpmethod', 'nearest', 'stdMethod', 'perColony');
% first round without normalization (to collect overall max and min of all), keep normalize = false
count = 0;
for posidx = posList
    condi = ceil(posidx/manMeta.posPerCondition);
    resOP = {};
    metaOP = meta;
    metaOP.conditions = meta.conditions(condi);
    metaOP.nWells = 1;
    metaOP.posPerCondition = 1;

    statsOP = cellStats(positions(posidx), metaOP, positions(1).dataChannels);

    clf
    options.conditionIdx = 1;
    options.colonyRadius = radii(condi);
    resOP = plotRadialProfiles(statsOP, metaOP, options);
    close all

    count = count+1;
    if count == 1
        maxvalsNuc = max(resOP.nuc_profile);
        minvalsNuc = min(resOP.nuc_profile);

        maxvalsCyt = max(resOP.cyt_profile);
        minvalsCyt = min(resOP.cyt_profile);

        maxvalsRatio = max(resOP.ratio_profile);
        minvalsRatio = min(resOP.ratio_profile);
    else
        maxvalsNuc = max(maxvalsNuc, max(resOP.nuc_profile));
        minvalsNuc = min(minvalsNuc, min(resOP.nuc_profile));

        maxvalsCyt = max(maxvalsCyt, max(resOP.cyt_profile));
        minvalsCyt = min(minvalsCyt, min(resOP.cyt_profile));

        maxvalsRatio = max(maxvalsRatio, max(resOP.ratio_profile));
        minvalsRatio = min(minvalsRatio, min(resOP.ratio_profile));
    end
end

if norMethodPos == 2
    nuclimits = [minvalsNuc' maxvalsNuc'];
    options.nuclimits = nuclimits;
    cytlimits = [minvalsCyt' maxvalsCyt'];
    options.cytlimits = cytlimits;
    ratiolimits = [minvalsRatio' maxvalsRatio'];
    options.ratiolimits = ratiolimits;
end

options.normalize = true; % set normalization here

for posidx = posList
    condi = ceil(posidx/manMeta.posPerCondition);
    options.conditionIdx = 1;
    options.colonyRadius = radii(condi);
    metaOP = meta;
    metaOP.conditions = meta.conditions(condi);
    metaOP.nWells = 1;
    metaOP.posPerCondition = 1;
    statsOP = cellStats(positions(posidx), metaOP, positions(1).dataChannels);
    figure()
    resOP = plotRadialProfiles(statsOP, metaOP, options);
    ylim([0 1.1]);
    xlim([0 radii(condi)]);
    axis square
    title([meta.conditions{condi} ' Pos' num2str(posidx)])
    saveas(gcf, ['radialProfileOP_' meta.channelLabel{:} '_' meta.conditions{condi} '_Pos' num2str(posidx) '.png'])
end

%% density plot condition by condition (MP Only)
rAll = cell(numel(meta.conditions),1);
densityAll = cell(numel(meta.conditions),1);

numberAll = zeros(manMeta.posPerCondition,numel(manMeta.conditions));
for condi = 1:numel(meta.conditions)
    rAll{condi} = res{condi}.r';
    density = [];
    for pos = 1:manMeta.posPerCondition
        posidx = (condi-1)*manMeta.posPerCondition+pos;
        densityTem = res{condi}.celldensity{posidx}';
        density = [density,densityTem];
        numberAll(pos,condi) = positions(posidx).ncells;
    end
    densityAll{condi} = density;
end

PI = 3.1415926;
densityColonyAll = numberAll/(PI*350^2)*10^6;
numberAllMean = mean(numberAll,1);
numberAllStd = std(numberAll,0,1);

c = lines(5);
fs = 20;
p = [];
figure()
for condi = 1:numel(meta.conditions)
    r = rAll{condi};
    densityMEAN = mean(densityAll{condi}*10^6,2); % cells/mm^2
    densitySTD = std(densityAll{condi}*10^6,0,2); % cells/mm^2
    hold on
    errorbar(r, densityMEAN, densitySTD,'--','LineWidth',0.5, 'Color', c(condi,:));
    ptemp = plot(r, densityMEAN,'-','LineWidth',3, 'Color', c(condi,:), 'DisplayName', manMeta.conditions{condi});
    hold off
    p = [p,ptemp];
end
xlim([0 350])
ylim([0 30000])
legend(p,'Location','best');
xlabel('edge distance ( um )')
ylabel('cell density (cells/mm^2)');
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')
set(gca, 'LineWidth', 2);
saveas(gcf, ['radialProfile_density' '.png'])

%% Calculate density & confluency of not confluent position
posNoConfluent = [1,5,9];
Area = zeros(numel(positions),1);
Confluency = ones(numel(positions),1);
Density = zeros(numel(positions),1);
for pi = posNoConfluent
    condi = ceil(pi / meta.posPerCondition);
    cpi = pi - (condi-1)*meta.posPerCondition;
    img = max(positions(pi).loadImage(dataDir, 0, 1),[],3);
    img = imadjust(img,stitchedlim(img));
    imgBW = imbinarize(img);
    imgBWCLEAN = bwareaopen(imgBW,1000); % remove small particles
    diskSize = 30; % play with disk size
    imgBW2 = imfill(imclose(imgBWCLEAN,strel('disk',diskSize)),'holes');
    figure()
    imshowpair(imgBWCLEAN,imgBW2,'montage')
    saveas(gcf,['Area_' manMeta.conditions{condi} '_Position_' num2str(pi) '.jpg']);
    Area(pi) = bwarea(imgBW2)*meta.xres*meta.yres/10^6;
end
for pi = 1:numel(positions)
    if Area(pi) == 0
        Area(pi) = stats.area*100;
    else
        Confluency(pi) = Area(pi)/(stats.area*100);
    end
    Density(pi) = positions(pi).ncells/Area(pi); % cells/mm^2
end
DensityAll = zeros(numel(manMeta.conditions),1);
DensitySTD = zeros(numel(manMeta.conditions),1);
coloniesPerSize = cellfun(@numel,filelist);
posbins = [0 cumsum(coloniesPerSize)];
for condi = 1:numel(manMeta.conditions)
    DensityAll(condi) = mean(Density(posbins(condi)+1:posbins(condi+1)));
    DensitySTD(condi) = std(Density(posbins(condi)+1:posbins(condi+1)));
end

%% scatter plot subpopulation on image

condi = 3; % condition index
pi = 2; % position index

cpi = pi + (condi-1)*meta.posPerCondition;
subDir = filelist{condi}(pi).folder;

channelThresholds = [0, stats.thresholds(2), stats.thresholds(3), stats.thresholds(4)]; % use automatic thresholds
%channelThresholds = [0 500 500 400]; % manually adjust thresholds

for ci = 1:numel(manMeta.channelLabel)
    % s = strsplit(positions(pi).filename,{'_FusionStitcher','.ims'});
    % prefix = [s{:}];
    %
    % zdir = [prefix '_zslices'];
    % img = imread(fullfile(subDir, zdir, sprintf([prefix '_MIP_w%.4d.tif'], ci-1)));
    img = max(positions(cpi).loadImage(dataDir,ci-1,1),[],3);
    % img = imadjust(img,stitchedlim(img));
    img = mat2gray(img, [150 4000]);
    [X,Y] = meshgrid(1:size(img,2),1:size(img,1));
    if strcmp(imageType,'MP')
        R = sqrt((X - positions(cpi).center(1)).^2 + (Y - positions(cpi).center(2)).^2);
        mask = R > positions(cpi).radiusPixel*1.05;
        img(mask) = max(img(:));
    end

    nucLevel = positions(cpi).cellData.nucLevel;
    background = positions(cpi).cellData.background;
    Ncells = positions(cpi).ncells;

    positive = {};
    for i = 1:meta.nChannels
        positive{i} = nucLevel(:,i) - background(i) > channelThresholds(i);
    end

    figure,
    imshow(img)
    hold on
    %scatter(positions(pi).cellData.XY(:,1),positions(pi).cellData.XY(:,2),'filled')
    XY = positions(cpi).cellData.XY(positive{ci},:);
    scatter(XY(:,1),XY(:,2),5,'filled','r')
    %XY = positions(pi).cellData.XY(AP2Cp,:);
    %scatter(XY(:,1),XY(:,2),'x','g')
    XY = positions(cpi).cellData.XY(positive{2},:);
    %scatter(XY(:,1),XY(:,2),'x','b')
    if strcmp(imageType,'MP')
        scatter(positions(cpi).center(1),positions(cpi).center(2),500,'.','g')
    end
    hold off
    title(meta.channelLabel{ci})
    saveas(gcf,['subpopulation_' meta.channelLabel{ci} '_Pos' num2str(cpi) '.jpg']);
end
% save new manually adjusted thresholds
for cpi = 1:numel(positions)
    positions(cpi).cellData.channelThresholds = channelThresholds;
end
save(fullfile(dataDir,'positions'), 'positions');
stats.thresholds = channelThresholds;
save(statsfile, 'stats');
%% radial probability of being positive condition by condition (MP Only)

Pall = {};
Pstdall = {};
xall = {};
for condi = 1:numel(meta.conditions)
    combos = {}; % [2 4],[2 3],[3 4]};
    stats.markerChannels = 2:4;
    [P,x,Pstd] = radialPositive(stats, positions, condi, meta, combos);
    Pall{condi} = P;
    Pstdall{condi} = Pstd;
    xall{condi} = x;
    ylim([0 1]);
    saveas(gcf, fullfile(dataDir, ['radialpositive_' meta.conditions{condi} '_', meta.channelLabel{2:end} '.png']));
end

%% radial probability of being positive channel by channel (MP Only)
markerChannels = 2:4;
for channel = markerChannels
    figure()
    colors = lines(numel(manMeta.conditions));

    lw = 2;
    fs = 20;
    legendstr = [];
    p = [];
    for i = 1:numel(meta.conditions)
        Pspe = Pall{i}(channel,:);
        Pstd_spe = Pstdall{i}(channel,:);
        xspe = xall{i};

        ptemp = plot(xspe,Pspe,'LineWidth',lw,'Color',colors(i,:));
        legendstr = [legendstr,{[meta.conditions{i},' ',meta.channelLabel{channel}]}];
        p = [p,ptemp];
        hold on
        good = ~isnan(Pspe);
        fill([xspe(good),fliplr(xspe(good))],...
            [Pspe(good) + Pstd_spe(good), fliplr(Pspe(good) - Pstd_spe(good))],...
            colors(i,:),'FaceAlpha',0.2,'EdgeColor','none');
    end
    xlim([0 350]);
    %ylim([0 1]);
    xlabel('edge distance ( um )')
    ylabel('positive fraction ');
    legend(p,legendstr,'FontSize',15,'Location','best')
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
    set(gca, 'LineWidth', 2);
    saveas(gcf, fullfile(dataDir, ['radialpositive_combine_', meta.channelLabel{channel} '.png']));
end
%% count populations

combo = [4 3 2];
conditionsidx = 1:numel(manMeta.conditions);
counts = countPopulations(positions, meta, stats, dataDir, combo, conditionsidx);
save(countfile, 'counts');

%% make pretty scatter plot

for condi = 1:numel(meta.conditions)

    close all;
    options = struct();

    options.selectDisType = [1, 1, 1, 1]; % 1 - nuc, 2 - cyto, 3 - NCratio, 4 - CNratio
    options.channelThresholds = stats.thresholds; % re-define channelThres when plot cyto or NCratio
    %options.checkPos = badidxCon; % whether you want to check defined cells in scatter plot

    options.imageType = imageType;
    options.conditionIdx = condi;
    options.channelCombos = {[4 3 2]}; % {[3 2 4], [4 3 2], [4 2 3]}

    options.axisLabel = {'DAPI','BRDU','BRA','YAP'};
    options.channelMax = exp([3 3 3 3])-1; %xlim or ylim
    options.log1p = [true true true true];
    if strcmp(imageType,'MP')
        options.radiusMicron = positions(posbins(condi)+1).radiusMicron;
    end
    options.conditionsCombined = false;
    %options.showThreshold = [false false true false];

    scatterMicropattern(stats, meta, dataDir, options)
end


%% pretty channel combo pie (MP Only)

options = struct();
options.channels = [2 4 3];
options.pieOrder = [3 1 2];
% tolerances in original order of channels
options.tol = [0.01 0.99; 0.01 0.99; 0.01 0.99; 0.01 0.99];
options.margin = 50;
options.scalebar = false;

posidx = [1 5 9]; % the first position sets the lookup table

for pi = posidx

    ci = ceil(pi/meta.posPerCondition);
    cpi = mod(pi-1, meta.posPerCondition)+1;
    fname = filelist{ci}(cpi).name;
    subDir = filelist{ci}(cpi).folder;

    if pi == posidx(1)
        Ilim = micropatternPieVis(subDir, positions(pi), options);
    else
        options.Ilim = Ilim;
        micropatternPieVis(subDir, positions(pi), options);
    end
end

%% multi condition pie (MP Only)

options = struct();
options.channels = [2 3 4];
% tolerances in original order of channels
options.tol = [0.01 0.99; 0.1 0.99; 0.01 0.99; 0.01 0.99];

conditionIdx = 1:numel(manMeta.conditions);

options.positionIdx = (conditionIdx-1)*meta.posPerCondition + 1;

micropatternPieVisConditions(dataDir, positions, options, meta);

