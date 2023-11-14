clear; close all; clc

%% load data, setup
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
dataDir = scriptPath;
savedir = fullfile(dataDir,'figures');
if ~exist(savedir,'dir'), mkdir(savedir); end

baseDir = fullfile(dataDir,'processed_data');
liveDir = fullfile(baseDir,'live');
fixedDir = fullfile(baseDir,'fixed');

%load metadata
liveMeta = load(fullfile(liveDir,'meta.mat'));
liveMeta = liveMeta.meta;

fixedMeta = load(fullfile(fixedDir,'combined_meta.mat'));
fixedMeta = fixedMeta.meta;
channelLabels = fixedMeta.channelLabel;

%load combined live and fixed cell quantification
load(fullfile(baseDir,'lt.mat'))
npos = length(lt.live_position);

%number of time points from live imaging
ntime = lt.live_position(1).nTime;
nucChannel = 0;
%vector of time points based on time interval from live imaging
treatmentTime = 3;
s = strsplit(liveMeta.timeInterval,'min');
tvec = ((1:ntime) - treatmentTime)*str2double(s{1})/60;

ppc = liveMeta.posPerCondition; %positions per condition
ncond = fixedMeta.nWells; %number of conditions
Rmax = lt.fixed_position(1).radiusMicron; %colony radius in microns

%labels for conditions for saving figures:
condlabels = strrep(liveMeta.conditions,',','_');
condlabels = strrep(condlabels,'+','_');

%plotting options:
lw = 3; %line width
fs = 28; %font size
lfs = 20; %font size for figure legends
figpos = figurePosition(560, 560); %figure position
reps = 100;

%% get fixed data in radial bins for each colony in each condition
%use the same number of radial bins for live and fixed data
nbin = 30;

%get channels with stains for cell fate markers
fixedchan = [2,4,7,8];
fixedfield = {'nucLevel','nucLevel','nucLevel','nucLevel'};
cls = fixedMeta.channelLabel(fixedchan); %channel labels

fixedopts = struct('fields',{fixedfield},'normalizeToDAPI',false(1,length(fixedchan)),...
    'nbin',nbin,'avgmethod','median');

Ys = cell(1,ncond); Rf = cell(1,ncond);
Yall = cell(1,npos); Rfall = cell(1,npos); Rall = cell(1,npos);
for cidx = 1:ncond
    pidxs = (cidx - 1)*ppc + (1:ppc);
    
    fixeddata = cell(ppc,1); rs = cell(ppc,1);
    for ii = 1:ppc
        pidx = pidxs(ii);
        [Y,r,~,edges] = newRadialProfile(lt.fixed_position(pidx),fixedchan,fixedopts);
        
        fixeddata{ii} = Y([1:end,end],:); rs{ii} = edges;
        Yall{pidx} = fixeddata{ii}; Rfall{pidx} = edges; Rall{pidx} = r;
    end
    Ys{cidx} = cell2mat(fixeddata);% Rf{cidx} = repmat(rs{1},ppc,1);
    Rf{cidx} = cell2mat(rs);
end

maxvals = cellfun(@(x) max(x,[],1),Ys,'UniformOutput',false);
maxvals = max(cell2mat(maxvals'),[],1);

minvals = cellfun(@(x) min(x,[],1),Ys,'UniformOutput',false);
minvals = min(cell2mat(minvals'),[],1);

Ynorm = cellfun(@(x) (x - minvals)./(maxvals - minvals),Ys,'UniformOutput',false);
Yall = cellfun(@(x) (x - minvals)./(maxvals - minvals),Yall,'UniformOutput',false);


%% plot overall profiles for each gene in each condition as fate maps
cmode = 'rgb';
chan = [1 2 3];
nc = length(chan);

%for each radial bin in each colony
Is = cell(1,npos);
for ii = 1:npos
    Y = Yall{ii}(:,chan);
    [~,I] = max(Y,[],2);
    Is{ii} = I;
    
end

profs = cell(1,ncond);
for cidx = 1:ncond
    prof = NaN(nbin+1,nc);
    pidxs = (cidx - 1)*ppc + 1:cidx*ppc;
    I = cell2mat(Is(pidxs));
    for ci = 1:nc
        prof(:,ci) = sum(I == ci,2)/ppc;
    end
    profs{cidx} = prof;
end

%patterns from radial profiles
opts = struct('Rmax',350,...
    'channels',{cls(chan)},...
    'order',[1 3 2],...
    'colorMode','rgb',...
    'rmode', 'edges');

for cidx = 1:2
    pidxs = (cidx - 1)*ppc + 1:cidx*ppc;
    edges = mean(cell2mat(Rfall(pidxs)),2); %Rfall{pidxs(1)}
%     [clusterim, clusterlabel, savelabel] =...
%         profile2pattern(profs{cidx},Rc{cidx},opts);
    [clusterim, clusterlabel, savelabel] =...
        profile2pattern(profs{cidx},edges,opts);
    figure
    imshow(clusterim{1})
    xlabel(clusterlabel)
    title(condlabels{cidx},'Interpreter','none')
    cleanSubplot
    drawnow
    
    imwrite(clusterim{1},fullfile(savedir,['fateMap_',condlabels{cidx},savelabel,'.png']))
end

%% plot aggregate signaling kymograph by condition
close all
livefield = {'NCratio'};
livechan = 2;
liveopts = struct('fields',{livefield},'normalizeToDAPI',false,'nbin',nbin,...
    'maxval',2,'minval',0,'avgmethod','median','suppressoutput',true);

X = NaN(nbin+1,ntime,ncond); R = NaN(nbin+1,ntime,ncond);
Err = NaN(nbin+1,ntime,ncond); T = repmat(tvec,nbin+1,1);
colors = turbo(nbin);
colors = colors(end:-1:1,:);
scale = 0.1;
for cidx = 1:ncond
    pidxs = (cidx - 1)*ppc + 1:cidx*ppc;
    condlabel = condlabels{cidx};
    disp(condlabel)
    for ti = 1:ntime
        liveopts.t = ti;
        [profiles,err,rad,edges] = newRadialProfile(lt.live_position(pidxs),livechan,liveopts);
        X(:,ti,cidx) = profiles([1:end,end]); R(:,ti,cidx) = edges; Err(:,ti,cidx) = err([1:end,end]);
    end
end


s4min = min(X(:,tvec>0,1),[],'all');
s4max = mean(X(1,tvec>0 & tvec > 9,1),'all');
X = (X - s4min)/(s4max - s4min);

for cidx = 1:2%ncond
    condlabel = condlabels{cidx};
    savename = strcat('kymographs_SMAD4_',condlabel,'.png');
    figure('Position',figpos)
    hold on
    surf(T,R(:,:,cidx),X(:,:,cidx),'LineStyle','none')
    cleanSubplot(fs); view(2); axis square; colormap turbo
    xlim([min(T(:)) max(T(:))]); ylim([min(R(:)) max(R(:))]);
    caxis([0.2,0.95])%caxis([0.4,0.9])
    xlabel('time (hr)'); ylabel('edge distance (um)')
    h = colorbar; h.Label.String = 'SMAD4 (N:C)';
    h.Location = 'northoutside'; h.TickLabels = '';
    
    drawnow
    savefigure(fullfile(dataDir,'figures',savename));
    
    %radial time series
    savename = strcat('timeseries_SMAD4_',condlabel,'.png');
    figure('Position',figpos); hold on
    for jj = 1:nbin
        plot(tvec,X(jj,:,cidx),'Color',colors(end-jj+1,:))
    end
    hold off
    xlim([min(tvec) max(tvec)]); ylim([-0.1,1.25])%ylim([0,1.3])%ylim([0.35,0.95])
    cleanSubplot(fs); axis square
    xlabel('time (hr)'); ylabel('SMAD4 (N:C)')
    drawnow
    savefigure(fullfile(savedir,savename));
end

%% signaling histories in radial bins in individual colonies
%previous block averaged colonies in the same condition, this block
%collects a set of radially averaged signaling histories for each
%individual colony

tic
Xs = cell(1,ncond); Rs = cell(1,ncond); Rt = cell(1,ncond);
Xall = cell(1,npos); Rall = cell(1,npos);
for cidx = 1:ncond
    pidxs = (cidx - 1)*ppc + (1:ppc);
    liveData = cell(ppc,1); rs = cell(ppc,1);
    for ii = 1:ppc
        pidx = pidxs(ii);
%         Xlive = NaN(nbin,ntime); R = NaN(nbin,ntime);
        Xlive = NaN(nbin+1,ntime); R = NaN(nbin+1,ntime);
        for ti = 1:ntime
            liveopts.t = ti;
            [X,~,rad,edges] = newRadialProfile(lt.live_position(pidx),livechan,liveopts);
%             Xlive(:,ti) = X; R(:,ti) = rad;
            Xlive(:,ti) = X([1:end,end]); R(:,ti) = edges;
        end
        Xlive = (Xlive - s4min)/(s4max - s4min);
        Rall{pidx} = R; Xall{pidx} = Xlive;
        liveData{ii} = NaN(size(Xlive,1),ntime); rs{ii} = NaN(size(Xlive,1),1);
        for bi = 1:size(Xlive,1)
            v = Xlive(bi,:); nanmask = isnan(v) | v > 2 | v < 0;
            if sum(~nanmask) > 2
                liveData{ii}(bi,:) =...
                    interp1(tvec(~nanmask),v(~nanmask),tvec,'linear','extrap');
                rs{ii}(bi) = R(bi,end);
            end
        end
    end

    Xs{cidx} = cell2mat(liveData); Rs{cidx} = repmat(rs{1},ppc,1);
end
toc

%% plot time series colored by distance from the colony edge
cidx = 1; %condition index: bmp4 without wnti
X = cell2mat(Xs(cidx)'); R = cell2mat(Rs(cidx)');

[~,I] = sort(R,'ascend');
nhists = length(I);

colors = turbo(nhists);
colors = colors(end:-1:1,:);

figure('Position',figurePosition(560,560))
hold on
for jj = 1:nhists
    plot(tvec,X(I(jj),:),'Color',colors(end-jj+1,:))
end
hold off
xlim([min(tvec) max(tvec)]);
ylim([0,1.4])
cleanSubplot(fs); axis square
xlabel('time (hr)')
ylabel('SMAD4 (N:C)')
savefigure(fullfile(savedir,['historiesByDist_',condlabels{cidx}]))

%% soft c means (fuzzy clustering)
k = 3;
cidx = 2;
chans = [1 3 2];
chanstr = cell2mat(strcat('_',cls(chans)));
order = [1 3 2];
colors = {'r','b','g'};

X = Xs{cidx}; R = Rs{cidx}; Y = Ynorm{cidx};
% [C, S, L,~,~,MU] = pca(X);
[~, S, ~,~,~,~] = pca(X);

ydata = Y(:,chans);
[centers,U] = fcm(X,k);
[~,Idx] = max(U,[],1);

[~, profs, rs, I] = reorderIndices(Idx',R,k);
U = U(I,:);
U = U(order,:);

opts = struct('Rmax',350,...
    'channels',{cellstr(strcat("cluster ",num2str((1:k)')))'},...
    'order',[1 3 2],...
    'colorMode','rgb',...
    'rmode', 'edges');
[clusterim, ~, ~] = profile2pattern(profs,rs,opts);

[profs, rs] = profsFromFuzzyClusters(U,R,k);
opts = struct('Rmax',350,...
    'channels',{cellstr(strcat("cluster ",num2str((1:k)')))'},...
    'order',[1 2 3],...
    'colorMode','rgb',...
    'rmode', 'edges');
[fuzzyim, ~, ~] = profile2pattern(profs,rs,opts);

figure('Position',figurePosition(2*560, 560))
subplot(1,2,1)
imshow(clusterim{1})
title('discrete assignments')
cleanSubplot(fs)

subplot(1,2,2)
imshow(fuzzyim{1})
title('fuzzy assignments')
cleanSubplot(fs)

if k == 3
    col = U';
elseif k == 2
    col = [U',zeros(size(U,2),1)];
end

imwrite(clusterim{1},fullfile(savedir,['cmeansMapDiscrete_',condlabels{cidx},sprintf('_k%d',k),'.png']))
imwrite(fuzzyim{1},fullfile(savedir,['cmeansMapFuzzy_',condlabels{cidx},sprintf('_k%d',k),'.png']))

%pca plot colored for cluster assignment
dim1 = 1; dim2 = 2;
figure('Position',figpos)
scatter(S(:,1),S(:,2),30,col,'filled')
cleanSubplot(fs); axis square
xlabel(sprintf('PC%d',dim1)); ylabel(sprintf('PC%d',dim2))
% title('color = clusters')
savefigure(fullfile(savedir,...
        ['PCA_cmeans_',sprintf('k%d_',k),condlabels{cidx},sprintf('_PC%d',[dim1 dim2])]))

%pca plot colored for fate
figure('Position',figpos)
scatter(S(:,1),S(:,2),30,ydata,'filled')
cleanSubplot(fs); axis square
xlabel(sprintf('PC%d',dim1)); ylabel(sprintf('PC%d',dim2))
% title('color = fate markers')
savefigure(fullfile(savedir,...
        ['PCA_fate',chanstr,'_',condlabels{cidx},sprintf('_PC%d',[dim1 dim2])]))

%plot cluster centers:
figure('Position',figpos); hold on
for ki = 1:k
    plot(tvec,centers(I(ki),:),'Color',colors{ki},'LineWidth',lw)
end
hold off
cleanSubplot(fs); axis square
xlabel('time (hr)')
ylabel('SMAD4 (N:C)')
xlim(tvec([1,end])); ylim([0.2,1])
savefigure(fullfile(savedir,['cmeansCentroids_',sprintf('k%d_',k),condlabels{cidx}]))

%plot histories colored by cluster membership:
figure('Position',figpos); hold on
for ii = 1:size(X,1)
    plot(tvec,X(ii,:),'Color',col(ii,:))
end
hold off
cleanSubplot(fs); axis square
xlabel('time (hr)')
ylabel('SMAD4 (N:C)')
xlim(tvec([1,end])); ylim([0,1.4])%ylim([0.5,1])
savefigure(fullfile(savedir,['cmeansAllHistories_',sprintf('k%d_',k),condlabels{cidx}]))

%% show results with different numbers of clusters
close all
%pick a condition, get live and fixed data, do PCA
cidx = 2;
condlabel = condlabels{cidx};
X = Xs{cidx}; R = Rs{cidx};
[~, S, ~,~,~,~] = pca(X);

%clustering
ks = 2:5; %numbers of clusters to test
figure('Position',figurePosition(560,560))
for ii = 1:length(ks)
    k = ks(ii);
    
    [~, U] = fcm(X,k);
    [~,Idx] = max(U,[],1);
    [idx, ~, ~, ~] = reorderIndices(Idx',R,k);
    subplot(2,2,ii); hold on
    for ki = 1:k
        scatter(S(idx == ki,1),S(idx == ki,2),'filled')
    end
    cleanSubplot(fs)
    xticklabels({}); yticklabels({})
    if ismember(ii,[1 3])
        ylabel('PC2')
    end
    if ismember(ii,[3 4])
        xlabel('PC1')
    end
    drawnow
    
end
savename = ['PCA_',condlabel,sprintf('_k%d',ks),'.png'];
savefigure(fullfile(savedir,savename))

ks = 2:5;
figure('Position',figurePosition(560,560))
clusterims = cell(length(ks),1);
for ii = 1:length(ks)
    k = ks(ii);
    [~, U] = fcm(X,k);
    [~,Idx] = max(U,[],1);
    [~, profs, rs, ~] = reorderIndices(Idx',R,k);
    [~,dprofs] = max(profs,[],2);

    opts = struct('Rmax',350,...
        'channels',{cellstr(strcat("cluster ",num2str((1:k)')))'},...
        'colors',lines(k),...
        'rmode', 'edges');
    
    [clusterim, ~] = profile2pattern_discrete(dprofs,rs,opts);
    clusterims{ii} = clusterim{1};
    subplot_tight(2,2,ii)
    imshow(clusterim{1})
    cleanSubplot
    drawnow
end

% combine cluster maps into a single pie plot
mhalf = round(size(clusterims{1},1)/2);
[X,Y] = meshgrid(1:size(clusterims{1},2),1:size(clusterims{1},1));
F = atan2((X - mhalf),-(Y - mhalf));

N = length(clusterims);

mask = cell(1,N); 
for ii = 1:N
    mask{ii} = ~imdilate(F < -pi + (ii-1)*2*pi/N | F > -pi + ii*2*pi/N,strel('disk',10));
end

lines = sum(cat(3,mask{:}),3) == 0;

As = zeros(size(clusterims{1}),class(clusterims{1}));
order = [3 1 2 4];
for zi = 1:size(clusterims{1},3)
    A = As(:,:,zi);
    for ii = 1:N
        cim = clusterims{order(ii)}(:,:,zi);
        A(mask{ii}) = cim(mask{ii});
    end
    A(lines) = 1;
    As(:,:,zi) = A;
end

figure; imshow(As)
savename = ['pieslice_',condlabel,sprintf('_k%d',ks),'.png'];
imwrite(As,fullfile(savedir,savename))

%% c means + IWP2
k = 2;
cidx = 1;
chans = [1 3 2];
chanstr = cell2mat(strcat('_',cls(chans)));
order = [1 2];
colors = {'r','g','b'};

X = Xs{cidx}; R = Rs{cidx}; Y = Ynorm{cidx};
% [C, S, L,~,~,MU] = pca(X);
[~, S, ~,~,~,~] = pca(X);
% S = (X-MU)*C;

ydata = Y(:,chans);
[centers,U] = fcm(X,k);
[~,Idx] = max(U,[],1);

[~, profs, rs, I] = reorderIndices(Idx',R,k);
U = U(I,:);
U = U(order,:);

if k == 3
    col = U';
elseif k == 2
    col = [U',zeros(size(U,2),1)];
end

opts = struct('Rmax',350,...
    'channels',{cellstr(strcat("cluster ",num2str((1:k)')))'},...
    'order',[1 2 0],...
    'colorMode','rgb',...
    'rmode', 'edges');
[clusterim, ~, ~] = profile2pattern(profs,rs,opts);

[profs, rs] = profsFromFuzzyClusters(U,R,k);
opts = struct('Rmax',350,...
    'channels',{cellstr(strcat("cluster ",num2str((1:k)')))'},...
    'order',[1 2 0],...
    'colorMode','rgb',...
    'rmode', 'edges');
[fuzzyim, ~, ~] = profile2pattern(profs,rs,opts);

figure('Position',figurePosition(2*560, 560))
subplot(1,2,1)
imshow(clusterim{1})
title('discrete assignments')
cleanSubplot(fs)

subplot(1,2,2)
imshow(fuzzyim{1})
title('fuzzy assignments')
cleanSubplot(fs)

imwrite(clusterim{1},fullfile(savedir,['cmeansMapDiscrete_',condlabels{cidx},sprintf('_k%d',k),'.png']))
imwrite(fuzzyim{1},fullfile(savedir,['cmeansMapFuzzy_',condlabels{cidx},sprintf('_k%d',k),'.png']))

%pca plot colored for cluster assignment
dim1 = 1; dim2 = 2;
figure('Position',figpos)
scatter(S(:,1),S(:,2),30,col,'filled')
cleanSubplot(fs); axis square
xlabel(sprintf('PC%d',dim1)); ylabel(sprintf('PC%d',dim2))
% title('color = clusters')
savefigure(fullfile(savedir,...
        ['PCA_cmeans_',sprintf('k%d_',k),condlabels{cidx},sprintf('_PC%d',[dim1 dim2])]))

%pca plot colored for fate
figure('Position',figpos)
scatter(S(:,1),S(:,2),30,ydata,'filled')
cleanSubplot(fs); axis square
xlabel(sprintf('PC%d',dim1)); ylabel(sprintf('PC%d',dim2))
% title('color = fate markers')
savefigure(fullfile(savedir,...
        ['PCA_fate',chanstr,'_',condlabels{cidx},sprintf('_PC%d',[dim1 dim2])]))

%plot cluster centers:
figure('Position',figpos); hold on
for ki = 1:k
    plot(tvec,centers(I(ki),:),'Color',colors{ki},'LineWidth',lw)
end
hold off
cleanSubplot(fs); axis square
xlabel('time (hr)')
ylabel('SMAD4 (N:C)')
xlim(tvec([1,end])); ylim([0,1.4])
savefigure(fullfile(savedir,['cmeansCentroids_',sprintf('k%d_',k),condlabels{cidx}]))

%plot histories colored by cluster membership:
figure('Position',figpos); hold on
for ii = 1:size(X,1)
    plot(tvec,X(ii,:),'Color',col(ii,:))
end
hold off
cleanSubplot(fs); axis square
xlabel('time (hr)')
ylabel('SMAD4 (N:C)')
xlim(tvec([1,end])); ylim([0,1.4])%ylim([0.5,1])
savefigure(fullfile(savedir,['cmeansAllHistories_',sprintf('k%d_',k),condlabels{cidx}]))

%% manually split up the branches in PCA space
cidx = 1;
X = cell2mat(Xs(cidx)'); R = cell2mat(Rs(cidx)'); Y = cell2mat(Ys(cidx)');
[~, S, ~,~,~,~] = pca(X);

k = 2;
xy = S(:,[1,2]);
figure; scatter(xy(:,1),xy(:,2),'filled')
xl = get(gca,'Xlim');
xm = mean(xl); xd = xl(2) - xl(1);
yl = get(gca,'Ylim');
ym = mean(yl); yd = yl(2) - yl(1);
xl = xm + 0.75*[-xd xd]; yl = ym + 0.75*[-yd yd];
xlim(xl); ylim(yl)

vs = cell(1,k);
roi = gobjects(k,1);
for ki = 1:k
    roi(ki) = drawpolygon;
end
title('adjust vertices')
pause

for ki = 1:k
    vs{ki} = roi(ki).Position;
end

savefigure(fullfile(savedir,['PCA_handdrawnboundaries_',condlabels{cidx}]))

figure('Position',figurePosition(560*k,560))
ins = cell(1,k);
for ki = 1:k
    ins{ki} = inpolygon(xy(:,1),xy(:,2),vs{ki}(:,1),vs{ki}(:,2));
    subplot(1,k,ki); hold on
    scatter(xy(:,1),xy(:,2),30,'k','filled')
    scatter(xy(ins{ki},1),xy(ins{ki},2),30,'r','filled')
    cleanSubplot; axis square
    xlabel('PC1'); ylabel('PC2')
end

savefigure(fullfile(savedir,['PCA_handdrawnmembership_',condlabels{cidx}]))

if k == 3
    Idxs = cell(1,4);
    test = NaN(size(X,1),1);
    test(ins{1}) = 1;
    test(ins{2} & ~ins{1} & ~ins{3}) = 2;
    test(ins{3}) = 3;
    Idxs{1} = test;

    test = NaN(size(X,1),1);
    test(ins{1}) = 1;
    test(ins{2} & ~ins{1}) = 2;
    test(ins{3} & ~ins{2}) = 3;
    Idxs{2} = test;

    test = NaN(size(X,1),1);
    test(ins{1} & ~ins{2}) = 1;
    test(ins{2} & ~ins{3}) = 2;
    test(ins{3}) = 3;
    Idxs{3} = test;

    test = NaN(size(X,1),1);
    test(ins{1} & ~ins{2}) = 1;
    test(ins{2}) = 2;
    test(ins{3} & ~ins{2}) = 3;
    Idxs{4} = test;
    
    order = [1 3 2];
elseif k == 2
    Idxs = cell(1,2);
    
    test = NaN(size(X,1),1);
    test(ins{1}) = 1;
    test(ins{2} & ~ins{1}) = 2;
    Idxs{1} = test;
    
    test = NaN(size(X,1),1);
    test(ins{1} & ~ins{2}) = 1;
    test(ins{2}) = 2;
    Idxs{2} = test;
    
    order = [1 2 0];
end

for ii = 1:size(Idxs,2)
    [idx, profs, rs, ~] = reorderIndices(Idxs{ii},R,k);
%     subplot(1,4,ii); hold on
    figure('Position',figpos); hold on
    for ki = 1:k
        mask = idx == ki;
        scatter(xy(mask,1),xy(mask,2),30,colors{ki},'filled')
    end
    cleanSubplot(fs); axis square
    xlabel('PC1'); ylabel('PC2')
    savefigure(fullfile(savedir,...
        ['PCA_handdrawn_',sprintf('k%d_',k),condlabels{cidx},sprintf('_split%d.png',ii)]))
    
    for ki = 1:k
        line(vs{ki}([1:end,1],1),vs{ki}([1:end,1],2),'LineWidth',2,'Color','k')
    end
    savefigure(fullfile(savedir,...
        ['PCA_handdrawn_',sprintf('k%d_',k),condlabels{cidx},sprintf('_split%d',ii),'_withbox.png']))
    
    opts = struct('Rmax',350,...
        'channels',{cellstr(strcat("cluster ",num2str((1:k)')))'},...
        'order',order,...
        'colorMode','rgb',...
        'rmode', 'radius');
    [handim, ~, ~] = profile2pattern(profs,rs,opts);
    imwrite(handim{1},fullfile(savedir,...
        ['fateHandDrawn_',sprintf('k%d_',k),condlabels{cidx},sprintf('_split%d.png',ii)]))
end

for ii = 1:size(Idxs,2)
    [idx, ~, ~, ~] = reorderIndices(Idxs{ii},R,k);
%     subplot(1,4,ii); hold on
    figure('Position',figurePosition(560,round(560/2*0.95))); hold on
    for ki = 1:k
        mask = idx == ki;
        scatter(xy(mask,1),xy(mask,2),30,colors{ki},'filled')
    end
    cleanSubplot(fs);% axis square
    xlabel('PC1'); ylabel('PC2')
    savefigure(fullfile(savedir,...
        ['PCA_handdrawn_',sprintf('k%d_',k),condlabels{cidx},sprintf('_split%d_narrow.png',ii)]))
    
    for ki = 1:k
        line(vs{ki}([1:end,1],1),vs{ki}([1:end,1],2),'LineWidth',2,'Color','k')
    end
    savefigure(fullfile(savedir,...
        ['PCA_handdrawn_',sprintf('k%d_',k),condlabels{cidx},sprintf('_split%d',ii),'_narrow_withbox.png']))
end

%% local functions

function [idx, profs, rs, I] = reorderIndices(Idx,R,k)
%organize the cluster indices in order of increasing average radial
%position
means = NaN(1,k);
for ki = 1:k
    means(ki) = mean(R(Idx == ki));
end
idx = NaN(size(Idx));
[~,I] = sort(means,'ascend');
for ki = 1:k
    idx(Idx == I(ki)) = ki;
end

rs = unique(R);
profs = zeros(length(rs),k);
for ii = 1:length(rs)
    val = sum(R == rs(ii));
    for ki = 1:k
        profs(ii,ki) =  sum(R == rs(ii) & idx == ki)/val;
    end
end

end

function [profs, rs] = profsFromFuzzyClusters(U,R,k)
%make radial profiles using membership values for each radial bin

if size(U,1) == k
    U = U';
end

%rescale for brighter colors
maxvals = max(U,[],1); minvals = min(U,[],1);
U = (U - minvals)./(maxvals - minvals);

rs = unique(R);
profs = zeros(length(rs),k);
for ii = 1:length(rs)
    for ki = 1:k
        profs(ii,ki) = mean(U((R == rs(ii)),ki));
    end
end

end






