clear; close all; 

scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);

% raw data location, modify if not the same as the location of this script
dataDir = scriptPath; 
cd(dataDir);

%manMeta = struct();
manMeta.nucChannel = 0;
manMeta.channelLabel = {'DAPI','SMAD2','PRDM1','pSMAD1'}; 

manMeta.conditions = {'100um','150um','200um','300um','700um'};
radii = [100 150 200 300 700]/2;

%manMeta.conditions = {'700um'};
manMeta.posPerCondition = 5; %minimum num colonies or files with the same condition
manMeta.nWells = numel(manMeta.conditions);
manMeta.nPositions = manMeta.posPerCondition*manMeta.nWells;

filelist = {};
filelist{1} = dir(fullfile(dataDir,'*PSMAD1_MP42h_100um*.ims'));
filelist{2} = dir(fullfile(dataDir,'*PSMAD1_MP42h_150um*.ims'));
filelist{3} = dir(fullfile(dataDir,'*PSMAD1_MP42h_200um*.ims'));
filelist{4} = dir(fullfile(dataDir,'*PSMAD1_MP42h_300um*Fusion*.ims'));
filelist{5} = dir(fullfile(dataDir,'*PSMAD1_MP42h_700um*Fusion*.ims'));

coloniesPerSize = cellfun(@numel,filelist);
posbins = [0 cumsum(coloniesPerSize)];
exclude = [49]; % bad colony - bright junk in PRDM1 channel

posfile = fullfile(dataDir,['positions' [manMeta.channelLabel{:} manMeta.conditions{:} '.mat']]);
if exist(posfile)
	load(posfile, 'positions');
    disp(['loading ' posfile]);
end

statsfile = fullfile(dataDir,['stats' [manMeta.channelLabel{:} '.mat']]);
if exist(statsfile)
	load(statsfile, 'stats');
    disp(['loading ' statsfile]);
    thresholds = stats.thresholds; % because the per condition stuff below regenerates stats
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

%% extractData
opts = struct('cytoSize',4,'cytoMargin',2,'IoU',0.5,'Zmethod','new');

positions(manMeta.nPositions) = Colony();

% for all conditions (aka wells)
tic
for ci = 1:manMeta.nWells
    
    % for all positions belonging to that condition
    for cpi = 1:manMeta.posPerCondition
    
        pi = manMeta.posPerCondition*(ci-1) + cpi;
        
        fname = filelist{ci}(cpi).name;
        subdir = filelist{ci}(cpi).folder;
        
        meta = Metadata(fullfile(subdir,fname), manMeta);
        % -1 is only because fusion stitcher attaches a black z-slice at the end
        % may have to be modified later
        meta.nZslices = meta.nZslices - 1;
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
        positions(pi) = Colony(meta, id);
        positions(pi).setRadius(radii(ci), meta.xres);
        positions(pi).well = ci;
        positions(pi).dataChannels = 1:positions(pi).nChannels;
        
        positions(pi).cellData = extractDataMultiZ(subdir, fname, meta, opts);
        positions(pi).setCenter();
        positions(pi).makeRadialAvgSeg();
        positions(pi).ncells = size(positions(pi).cellData.XY,1);
    end
end
toc
save(posfile,'positions');

% %% combine positions colonies
% 
% pos100 = load( 'positionsDAPISMAD2PRDM1PSMAD1100um.mat', 'positions' );
% pos150 = load( 'positionsDAPISMAD2PRDM1PSMAD1150um.mat', 'positions' );
% pos200 = load( 'positionsDAPISMAD2PRDM1PSMAD1200um.mat', 'positions' );
% pos300 = load( 'positionsDAPISMAD2PRDM1PSMAD1300um.mat', 'positions' );
% pos700= load( 'positionsDAPISMAD2PRDM1PSMAD1700um.mat', 'positions' );
% 
% colonies = {pos100, pos150, pos200, pos300, pos700};
% 
% sumcol = 0;
% for i = 1:length(coloniesPerSize)
%     sumcol = sumcol + coloniesPerSize(i);
% end
% 
% totalpositions = pos100.positions;
% count=1; 
% 
% for i = 1:length(coloniesPerSize)
%   for j = 1:coloniesPerSize(i)
%         totalpositions(count) = colonies{i}.positions(j); 
%         count = count +1; 
%    end 
% end

%save(totalpos,'totalpositions');

%% make xsections of segmentation
%save 1 of these per condition
tic
for ci = 1:manMeta.nWells
    cpi = 1;
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
    nuclearim = squeeze(img(:,:,nucChannel+1,1:meta.nZslices));
    %change yinds and xinds if preferred
    yinds = 1400;
    xinds = 1:size(img,2);
    makeCrossSection(nuclearim,LM,meta,yinds,xinds)
    saveas(gcf, fullfile(dataDir, [prefix  'segmentationXsection_y' num2str(yinds) '.png'])); 
    close all
end
toc

%% stats for each size (which have different numbers of colonies)

stats_percond = {};

for condi = 1:numel(meta.conditions)
    
    metaX = meta;
    metaX.conditions = meta.conditions(condi);
    metaX.nWells = 1;
    posidx = posbins(condi)+1:posbins(condi+1);
    posidx = setdiff(posidx, exclude);
    metaX.posPerCondition = numel(posidx);
    
    stats_tmp = cellStats(positions(posidx), metaX, positions(1).dataChannels);
    stats_tmp.radiusMicron = {radii(condi)};
    stats_tmp.radiusPixel = {round(radii(condi)/meta.xres)};
    stats_tmp.thresholds = thresholds;
    
    % CLEAN UP
    % throwout brightest DAPI (dying/dead cells)
    i = 1;
    badidx = stats_tmp.nucLevel{1}(:,i) > mean(stats_tmp.nucLevel{1}(:,i)) + 2*std(stats_tmp.nucLevel{1}(:,i));
    % throw out darkest DAPI (usually segmentation error)
    badidx = badidx | (stats_tmp.nucLevel{1}(:,i) < mean(stats_tmp.nucLevel{1}(:,i)) - 2*std(stats_tmp.nucLevel{1}(:,i)) );
    % throw out cells outside the colony edges
    %badidx = badidx | sqrt(sum(stats.XY{condi}.^2,2))*meta.xres > radii(condi);
    
    stats_tmp.nucLevel{1} = stats_tmp.nucLevel{1}(~badidx,:);
    stats_tmp.cytLevel{1} = stats_tmp.cytLevel{1}(~badidx,:);
    stats_tmp.XY{1} = stats_tmp.XY{1}(~badidx,:);
    stats_tmp.sample{1} = stats_tmp.sample{1}(~badidx,:);
    
    if condi == 1
        stats = stats_tmp;
    else
        stats.conditions = [stats.conditions, stats_tmp.conditions];
        stats.radiusPixel = [stats.radiusPixel, stats_tmp.radiusPixel];
        stats.radiusMicron = [stats.radiusMicron, stats_tmp.radiusMicron];
        stats.nucLevel = [stats.nucLevel, stats_tmp.nucLevel];
        stats.cytLevel = [stats.cytLevel, stats_tmp.cytLevel];
        stats.XY = [stats.XY, stats_tmp.XY];
        stats.sample = [stats.sample, stats_tmp.sample];
        stats.othercelldata = [stats.othercelldata, stats_tmp.othercelldata];
    end
    
    stats_percond{condi} = stats_tmp;
end

tolerance = 0.01;
nbins = 50;
stats.makeHistograms(nbins, tolerance);
%save(statsfile, 'stats');

% clf
% conditionIdx = 1:numel(meta.conditions);
% for channelIdx = 1:4
%     options = struct('channelIndex',channelIdx, 'cytoplasmic', false,...
%                                 'cumulative',false, 'time', 1,...
%                             'conditionIdx',conditionIdx,...
%                             'titlestr',meta.channelLabel{channelIdx}...
%                             );
%     clf                
%     stats.plotDistributionComparison(options)
%     saveas(gcf, fullfile(dataDir, [meta.channelLabel{channelIdx(1)} '_dist.png'])); 
% end
% %ylim([0 0.02])
% %xlim([0 1000])

%% NCR for SMAD2
% not informative, we don't seem to get a good cytoplasmic stain

condi = 5;

R = sqrt(sum(stats.XY{condi}.^2,2))*meta.xres;
nL = stats.nucLevel{condi};
cL = stats.cytLevel{condi};
NCR = nL./cL;

scatter(R, NCR(:,2),'.')
ylim([0.5 3]);
%scatter(R, nL(:,2),'.r')

%% make intensity radial profile

% make intensity radial profile NEW

options = struct('nucChannels', 1:3,...
                'normalize', false, 'std', true, 'FontSize', 15, 'legend',true,...%);
                   'pointsPerBin', 400);
res = {};

% first round with normalization (to collect overall max and min)
figure,
for condi = 1:numel(meta.conditions)
    clf
    options.conditionIdx = condi;
    options.colonyRadius = radii(condi);
    res{condi} = plotRadialProfiles(stats, meta, options);
    if condi == 1
        maxvals = max(res{condi}.nuc_profile);
        minvals = min(res{condi}.nuc_profile);
    else
        maxvals = max(max(res{condi}.nuc_profile), maxvals);
        minvals = min(min(res{condi}.nuc_profile), minvals);
    end
end

% second round with normalization 
nuclimits = [minvals' maxvals'];
options.nuclimits = nuclimits;
options.normalize = true;

for condi = 1:numel(meta.conditions)
    clf;
    options.conditionIdx = condi;
    options.colonyRadius = radii(condi);
    res{condi} = plotRadialProfiles(stats, meta, options);
    ylim([0 1.1]);
    xlim([0 radii(condi)]);
    axis square
    title(meta.conditions{condi})
    saveas(gcf, ['radialProfile_' meta.channelLabel{:} '_' meta.conditions{condi} '.png']) 
end

%%
% combine different colony sizes in one plot

figure,
for ci = 2:4
clf
hold on
colors = lines(numel(meta.conditions));
for condi = 1:numel(meta.conditions)
    v(condi) = min(res{condi}.r);
    x = res{condi}.r;% - v(condi); %
    y = res{condi}.nuc_profile(:,ci);
    plot(x, y,'-', 'LineWidth',4, 'Color', colors(condi,:))
end
for condi = 1:numel(meta.conditions)
    x = res{condi}.r;% - min(res{condi}.r); %- v(condi);% 
    y = res{condi}.nuc_profile(:,ci)';
    e = res{condi}.nuc_profile_std(:,ci)';
    %e = res{condi}.nuc_profile_colstd(:,ci)';
    % shaded error bars
    fill([x, fliplr(x)],[(y + e), fliplr((y - e))], colors(condi,:),'FaceAlpha',0.1,'EdgeColor','none');
end
hold off
legend(meta.conditions)
title(meta.channelLabel{ci})
xlim([0 200]);
ylim([0 1.2]);

% make it pretty
fs = 32;
fgc = 'k';
bgc = 'w';
graphbgc = 1*[1 1 1]; 

xlabel('edge distance ( um )', 'FontSize',fs,'FontWeight','Bold','Color',fgc)
ylabel('intensity (a.u.)', 'FontSize',fs,'FontWeight','Bold','Color',fgc);

set(gcf,'color',bgc);
set(gca, 'LineWidth', 2);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')
set(gca,'XColor',fgc);
set(gca,'YColor',fgc);
set(gca,'Color',graphbgc);
axis square

saveas(gcf,['sizesCombined_' meta.channelLabel{:} '_' meta.channelLabel{ci} '.png'])
end


%% make intensity radial profile OLD

options = struct('nucChannels', 1:3,...
                'normalize', true, 'std', true, 'FontSize', 15, 'legend',true, ...
                    'nuclimits', nuclimits);
                
res = {};

options.individualColonies = false;
options.normalize=false;

for condi = 1:numel(meta.conditions)
    
    metaX = meta;
    metaX.posPerCondition = coloniesPerSize(condi);
    metaX.conditions = meta.conditions(condi);
    metaX.nWells = 1;
    posidx = posbins(condi)+1:posbins(condi+1);

    figure('Position',[0 0 500 500])
    res{condi} = plotRadialProfiles_old(positions(posidx), metaX, options);
    %ylim([0 1.1]);
    axis square
    title(meta.conditions{condi})
    saveas(gcf, ['radialProfileOld_' meta.channelLabel{:} '_' meta.conditions{condi} '.png']) 
end

%%

% combine different colony sizes in one plot
figure,
ci = 2;
clf
hold on
for condi = 1:numel(meta.conditions)
    errorbar(res{condi}.r, res{condi}.nuc_profile(:,ci), res{condi}.nuc_profile_std(:,ci),...
        'LineWidth',2)
end
hold off
legend(meta.conditions)
xlim([-20 150]);
ylim([0 1.2]);

%% thresholds

confidence = 0.95;
whichthreshold = [1 1 2];
conditions = 1:numel(meta.conditions);
stats.getThresholds(confidence, conditions);%, whichthreshold);

%save(statsfile, 'stats');

%% SAMPLE QC PLOT
% check that scatterplots of different samples in same conditions look
% similar

condi = 5;
i = 2;
j = 3;

X = stats.nucLevel{condi}(:,i);
Y = stats.nucLevel{condi}(:,j);
X = log(1+X/mean(X));
Y = log(1+Y/mean(Y));
scatter(X, Y, 2, stats.sample{condi})
xlim([0 1.5]);
ylim([0 3]);

%%
% check area
i = 1;
condi = 5;
hist(stats.othercelldata{condi}.nucArea,100)

% DONE:
% extend scatterMicropattern to more than 3 channels
% filter on stats
% include all cellData in stats
% make plotRadialProfiles return the profiles (to combine externally)

% discuss with Seth the problems with measuring area

%% scatter plot subpopulation on image

figure, 

normIdx = 1;
condi = 5;
cpi = 3; % [5,3] has strong artefact in psmad1 radial profile
pi = posbins(condi) + cpi; %position index
ci = 3; %channel index

% when specifying pi
%condi = ceil(pi / meta.posPerCondition);
%cpi = pi - (condi-1)*meta.posPerCondition;

channelThresholds = stats.thresholds;

% LOAD IMAGE
s = strsplit(positions(pi).filename,{'_FusionStitcher','.ims'});
prefix = [s{:}];
subDir = filelist{condi}(cpi).folder;
zdir = [prefix '_zslices'];
img = imread(fullfile(subDir, zdir, sprintf([prefix '_MIP_w%.4d.tif'], ci-1)));
%img = max(positions(pi).loadImage(dataDir,ci-1,1),[],3);
img = imadjust(img,stitchedlim(img));
%img = mat2gray(img, [150 4000]);

% MASK IMAGE
[X,Y] = meshgrid(1:size(img,2),1:size(img,1));
R = sqrt((X - positions(pi).center(1)).^2 + (Y - positions(pi).center(2)).^2);
mask = R > positions(pi).radiusPixel*1.05;
img(mask) = max(img(:));

% COLLECT CELL DATA
% from stats so including filters
idx = stats.sample{condi} == cpi-1;
XY = stats.XY{condi}(idx,:);
nucLevel = stats.nucLevel{condi}(idx,:);
area = stats.othercelldata{condi}.nucArea(idx,:);

% directly
%nucLevel = positions(pi).cellData.nucLevel;
%background = positions(pi).cellData.background;
%Ncells = positions(pi).ncells;

% SELECT CELLS TO VISUALIZE
positive = {};
for i = 1:meta.nChannels
    %positive{i} = nucLevel(:,i) - background(i) > 0.1*channelThresholds(i);
    positive{i} = nucLevel(:,i) > channelThresholds(i);
end
% posidx = (nucLevel(:,3)/channelThresholds(3) > 0.3) & (nucLevel(:,4)/channelThresholds(4) < 0.1);
% posidx = area < mean(area) - std(area); 
% posidx = positive{4} < channelThresholds(4);
% posidx = (nucLevel(:,3)/channelThresholds(3) < 1)...
%             & (nucLevel(:,4)/channelThresholds(4) < 0.5)...
%             & sqrt(sum(XY.^2,2))*meta.xres > 330;
posidx = nucLevel(:,1) > 0;

% figure,
% scatter(nucLevel(:,3)/channelThresholds(3),nucLevel(:,4)/channelThresholds(4),[],posidx)
% xlim([0 3]);
% ylim([0 1.5]);

% VISUALIZE
%figure,
imshow(img)
hold on
%scatter(positions(pi).cellData.XY(:,1),positions(pi).cellData.XY(:,2),'filled')
XYpos = XY(posidx,:);
% scatter(XYpos(:,1) + positions(pi).center(1),...
%         XYpos(:,2) + positions(pi).center(2),20,'filled','r')

% scatter(XYpos(:,1) + positions(pi).center(1),...
%         XYpos(:,2) + positions(pi).center(2),20,'filled','g')

myc = mat2gray(positive{ci});
%myc = imadjust(mat2gray(nucLevel(:,ci)));
%myc = imadjust(mat2gray(area(posidx)));
c = repmat(myc,1, 3);
c(:,3)=c(:,3)*0;
scatter(XYpos(:,1) + positions(pi).center(1),...
        XYpos(:,2) + positions(pi).center(2),30,c,'filled')

scatter(positions(pi).center(1),positions(pi).center(2),1000,'.','g')
hold off
title(meta.channelLabel{ci})

% for pi = 1:numel(positions)
%     positions(pi).cellData.channelThresholds = channelThresholds;
% end
% save(fullfile(dataDir,'positions'), 'positions');


%% SMAD1 levels do not go to zero on the colony edge
% so there was a problem with making radial profiles

ci = 2;
figure,
r = sqrt(sum(XYpos.^2,2))*meta.xres;
scatter(r,nucLevel(:,ci),'.')
hold on 
plot(positions(pi).radialProfile.BinEdges(1:end-1)*meta.xres,...
        positions(pi).radialProfile.NucAvgSeg(:,ci),'LineWidth',2,'Color','g')
hold off

% N = 50;
% radprof = zeros([N-1 1]);
% binranges = linspace(0,350,N);
% [bincounts,ind]= histc(r,binranges);
% for i = 1:N-1
%     radprof(i) = mean(nucLevel(ind==i,ci));
% end
% hold on
% plot(binranges(1:end-1), radprof,'LineWidth',2)
% j = 49;
% sum(ind==j)
% scatter(r(ind==j),nucLevel(ind==j,ci),'.','g')
% hold off

% average with equal number of datapoints per 'bin'
[B,I] = sort(r);

cellsperbin = 100;
N = round(numel(I)/cellsperbin);
edges = round(linspace(1, numel(I), N));

rprof = zeros([N-1 1]);
radprof = zeros([N-1 4]);

for i = 1:N-1
    rprof(i) = mean(r(I(edges(i):edges(i+1))));
    radprof(i,:) = mean(nucLevel(I(edges(i):edges(i+1)),:));
end
hold on
plot(rprof,radprof(:,ci),'LineWidth',2,'Color','r')
hold off

ylim([0 2000])
% % average with equal number of datapoints per 'bin' for all colonies
% % combined
% r = sqrt(sum(stats.XY{condi}.^2,2))*meta.xres;
% nL = stats.nucLevel{condi};
% 
% % average with equal number of datapoints per 'bin'
% [B,I] = sort(r);
% 
% ptsperbin = 200;
% N = round(numel(I)/ptsperbin);
% edges = round(linspace(1, numel(I), N));
% 
% rprof = zeros([N-1 1]);
% radprof = zeros([N-1 4]);
% 
% for i = 1:N-1
%     rprof(i) = mean(r(I(edges(i):edges(i+1))));
%     radprof(i,:) = mean(nL(I(edges(i):edges(i+1)),:));
% end
% hold on
% plot(rprof,radprof(:,ci),'LineWidth',2)

%% radial probability of being positive

Pall = {};
xall = {};
%channelThresholds = [0 500 500 800];
%stats.thresholds =  channelThresholds;%[0 500 400 700];
markerChannels = 2:4;

% NOTE : this is a mess, the first stats_percond is actually the full
% stats object but it works because we're indexing the first conditions of
% it, then all the other objects are actually

for condi = 1:numel(meta.conditions)

    combos = {};%[2 3],[2 4],[3 4]};
    metaX = meta;
    posidx = setdiff(posbins(condi)+1:posbins(condi+1), exclude);
    metaX.posPerCondition = numel(posidx);
    stats_percond{condi}.thresholds = stats.thresholds;
    stats_percond{condi}.markerChannels = 2:4;
    [P,x,Pstd] = radialPositive(stats_percond{condi}, positions(posidx), 1, metaX, combos);
    Pall{condi} = P;
    Pstdall{condi} = Pstd;
    xall{condi} = x;
    ylim([0 1]);
    %saveas(gcf, fullfile(dataDir, ['radialpositive_' meta.conditions{condi} '_', meta.channelLabel{2:end} '.png']));
    %close;
end

% problems: using the density estimate from all data instead of per colony
% getting the range right for every condition

%% combine different colony size in one plot

ci = 3;
figure,
hold on
for condi = 1:numel(meta.conditions)
    %plot(radii(condi) - xall{condi}*meta.xres, Pall{condi}(ci,:),'LineWidth',2)
    colRadius = radii(condi);
    xi = xall{condi};
    P = Pall{condi};
    x = colRadius - xi*meta.xres - v(condi) + v(5);
    plot(x, P(ci,:),'LineWidth',4,'Color',colors(condi,:))
end
for condi = 1:numel(meta.conditions)
    %plot(radii(condi) - xall{condi}*meta.xres, Pall{condi}(ci,:),'LineWidth',2)
    colRadius = radii(condi);
    xi = xall{condi};
    P = Pall{condi};
    Pstd = Pstdall{condi};
    good = ~isnan(P(ci,:));
    x = colRadius - meta.xres*xi(good) - v(condi) + v(5);
	fill([x,fliplr(x)],[P(ci,good) + Pstd(ci,good), fliplr(P(ci,good) - Pstd(ci,good))],colors(condi,:),'FaceAlpha',0.1,'EdgeColor','none');
end
legend(meta.conditions)
xlim([0 150]);
ylim([0 1]);
axis square
hold off
title(meta.channelLabel{ci})

% make it pretty
fs = 32;
fgc = 'k';
bgc = 'w';
graphbgc = 1*[1 1 1]; 

xlabel('edge distance ( um )', 'FontSize',fs,'FontWeight','Bold','Color',fgc)
ylabel('positive fraction', 'FontSize',fs,'FontWeight','Bold','Color',fgc);

set(gcf,'color',bgc);
set(gca, 'LineWidth', 2);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')
set(gca,'XColor',fgc);
set(gca,'YColor',fgc);
set(gca,'Color',graphbgc);


saveas(gcf, fullfile(dataDir, ['radialpositive_' meta.conditions{:} '_', meta.channelLabel{ci} '.png']));

%% count populations (doesn't work)
            
combo = [3 2 4];
%combo = [4 2];
conditionsidx = 1:5;%6:9;
combosubset = [1 3];
counts = countPopulations(positions, meta, stats, dataDir, combo, conditionsidx);
save(countfile, 'counts');

%% make pretty scatter plot

% since scatterplot is normalized in threshold units, for the smads set the
% threshold to the maximal value of the normalization of radialprofiles so
% they are in the same units
stats.thresholds(1) = maxvals(1);
stats.thresholds(2) = maxvals(2)*0.65;
stats.thresholds(4) = maxvals(4)*0.2;

for condi = [1 5]

close all;
options = struct();
options.conditionIdx = condi;%:numel(meta.conditions);
options.channelCombos = {[2 4 3]};% [2 4 3] 
options.channelThresholds = stats.thresholds;
options.norm = cat(3, [minvals(1) minvals(2) 0 minvals(4)],...
                    [maxvals(1) maxvals(2) stats.thresholds(3) maxvals(4)]);

options.axisLabel = {'DAPI','nuclear SMAD2/3', 'PRDM1', 'pSMAD1'};
options.channelMax = exp([3 1.5 2.5 1.25])-1; %xlim or ylim 
options.log1p = [true true true true];
options.radiusMicron = positions(posbins(condi)+1).radiusMicron;
options.conditionsCombined = false;
%options.showThreshold = [false false true false]; 

scatterMicropattern(stats, meta, dataDir, options)
end

%% see if fraction of cells with right SMAD1/2 range correlates with PRDM1+

x = [];
y = [];

for condi = 1:5
    nucLevel = stats.nucLevel{condi};

    % SELECT CELLS TO VISUALIZE
    positive = {};

    positiveS2 = nucLevel(:,2) > 0.5*stats.thresholds(2) & nucLevel(:,2) < 1.2*stats.thresholds(2);
    positiveS1 = nucLevel(:,4) > 0.5*stats.thresholds(4) & nucLevel(:,4) < stats.thresholds(4);
    x(condi) = sum(positiveS1 & positiveS2)/numel(positiveS1)
end

for condi = 1:5
    nucLevel = stats.nucLevel{condi};

    % SELECT CELLS TO VISUALIZE
    positive = {};

    positivePGC = nucLevel(:,3) > 1*stats.thresholds(3);
    y(condi) = sum(positivePGC)/numel(positivePGC)
end

plot(x,y,'-x')
xlabel('smad1/2 ok');
ylabel('prdm1+');


%% pretty channel combo pie
   
options = struct();
options.channels = [2 4 3];
options.pieOrder = [3 1 2];
% tolerances in original order of channels
% manMeta.channelLabel = {'DAPI','AP2C','NANOG','PRDM1'};
%manMeta.channelLabel = {'DAPI','AP2C','EOMES','SOX17'};
options.tol = [0.01 0.99; 0.3 0.99; 0.6 0.995; 0.1 0.99];% 0.3 0.99
%posidx = [6];% 5 7 8]; % the first position sets the lookup table
options.margin = 50;
options.scalebar = false;

for ci = 1
    for cpi = 15 % coloniesPerSize(ci)-1
    
    pi = posbins(ci) + cpi;

    %for pi = posidx

        %ci = ceil(pi/meta.posPerCondition);
        %cpi = mod(pi-1, meta.posPerCondition)+1;
        fname = filelist{ci}(cpi).name;
        subDir = filelist{ci}(cpi).folder;

        %if pi == posidx(1)
            Ilim = micropatternPieVis(subDir, positions(pi), options);
            Ilim{3}
        %else
            options.Ilim = Ilim;
            micropatternPieVis(subDir, positions(pi), options);
        %end
    %end
    end
end

%% multi condition pie
  
options = struct();
options.channels = [2 3 4];
% tolerances in original order of channels
% manMeta.channelLabel = {'DAPI','AP2C','NANOG','PRDM1'};
options.tol = [0.01 0.99; 0.1 0.99; 0.01 0.99; 0.7 0.995];
 
conditionIdx = [1 2 4 5 6];

options.positionIdx = (conditionIdx-1)*meta.posPerCondition + 1;

micropatternPieVisConditions(dataDir, positions, options, meta);



