clear; close all; clc

%% load data
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
dataDir = fullfile(scriptPath,'data');

T = struct();

%read sequencing count data
T.dose_cpm = readtable(fullfile(dataDir,'Full_Analysis_Data_7230.txt'));
T.time_cpm = readtable(fullfile(dataDir,'Full_Analysis_Data_7072.txt'));

T.dose_fpkm = readtable(fullfile(dataDir,'7230_gene_FPKM.annot.txt'));
T.time_fpkm = readtable(fullfile(dataDir,'7072_gene_FPKM.annot.txt'));

T.dose_fpkm.Symbol = T.dose_fpkm.external_gene_name;
T.time_fpkm.Symbol = T.dose_fpkm.external_gene_name;

fs = 28; %font size
lfs = 20; %legend font size
ms = 10; %marker size
figpos = figurePosition(560,560); %figure size
tol = 0.001; %tolerance
lw = 3;

%load signaling data
load(fullfile(dataDir,'signalData.mat'))
%average signaling level in each condition
SMAD4 = mean(X(6:end,:),1);

%list of transcription factor names and identifiers (not in the same order)
TFs = readcell(fullfile(dataDir,'TF_names_v_1.01.txt'));
TFEnsembl = readcell(fullfile(dataDir,'TFs_Ensembl_v_1.01.txt')); %<- use this one

%% parse data into usable matrices
% T.Properties.VariableNames

dsets = {'dose_cpm','time_cpm','dose_fpkm','time_fpkm'};

%just the values in each condition in a matrix
M = struct();
M.dose_cpm = table2array(T.dose_cpm(:,6:11));
M.time_cpm = table2array(T.time_cpm(:,[12 6:11]));
M.time_cpm = [mean(M.time_cpm(:,1:2),2),M.time_cpm(:,3:end)];
% M.dose_cpm_lin = exp(log(2)*M.dose_cpm);
% M.time_cpm_lin = exp(log(2)*M.time_cpm);
M.dose_cpm_lin = 2.^(M.dose_cpm);
M.time_cpm_lin = 2.^(M.time_cpm);

M.dose_fpkm_lin = table2array(T.dose_fpkm(:,5:10));
M.time_fpkm_lin = table2array(T.time_fpkm(:,[11 5:10]));
M.time_fpkm_lin = [mean(M.time_fpkm_lin(:,1:2),2),M.time_fpkm_lin(:,3:end)];
M.dose_fpkm = log2(0.1 + M.dose_fpkm_lin);
M.time_fpkm = log2(0.1 + M.time_fpkm_lin);

% times = [-0.5, 0, 4, 12, 20, 28, 42];
times = [0, 4, 12, 20, 28, 42];
doses = [0, 10, 30, 100, 300, 1000];

% CPM = log transform is log2(2 + reads) / (mean(library size)/10^6)
% for CPM : TMM normalization instead of just library size,
% adds weight to library size for adjustment between samples
% for FPKM no adjustment, but divide by transcript length as well

% FPKM = reads / (mean(library size)*(transcript length)/10^6)
% logFPKM = log2(0.1 + FPKM)

% 2^M = (E^log(2))^M = E^log(2)M


%% time series -> hierarchical clustering
close all
%methods: average, complete, single
%methods only appropriate for euclidean distance: centroid, ward, median
method = 'ward';

%dist: euclidean, spearman, correlation, cosine
dist = 'euclidean';

dset = 'time_cpm_lin';
Ydata = M.(dset);

%remove genes with names beginning in the character '.' (novel transcripts
%and pseudogenes)
if strcmp(dset,'time_cpm_lin') || strcmp(dset,'time_cpm')
    genes = T.time_cpm.Symbol;
    foldChanges = T.time_cpm.RowSums;
elseif strcmp(dset,'time_fpkm_lin') || strcmp(dset,'time_fpkm')
    genes = T.time_fpkm.Symbol;
    foldChanges = []; %<- need to calculate this for the fpkm stuff but not that hard
end

ydata = Ydata./max(Ydata,[],2); %normalize to maximum for plotting
% plot sets of genes of interest
IFgenes = {{'ISL1','TFAP2C','GATA3','HAND1'},...
    {'TFAP2A','GATA2','GABRP','KRT7','TP63','CDX2'},{'POU5F1','SOX2','NANOG'},{'ID1','ID2','ID3','ID4'}};
savelabels = {'amnion','amnion_2','pluri','IDs'};
legposes = [0.2315 0.67 0.277 0.248; 0.2327 0.5848 0.2768 0.3688; 0.6304 0.5589 0.2768 0.1875; 0.7494 0.2438 0.1661 0.2479];
for ii = 1:length(IFgenes)
    figure('Position',figurePosition(560,560)); hold on
    for jj = 1:length(IFgenes{ii})
        plot(times,ydata(strcmp(genes,IFgenes{ii}{jj}),:),'-x',...
            'LineWidth',lw,'MarkerSize',ms)
    end
    hold off
%     plot(times,ydata(ismember(genes,IFgenes{ii}),:),'LineWidth',2)
    cleanSubplot(fs); axis square
    xlabel('time (hr)'); ylabel('normalized CPM')
    legend(IFgenes{ii},'FontSize',lfs,'FontWeight','bold','Position',legposes(ii,:));
    ylim([0,1])
    savename = ['timeseriesexamples_',dset,'_',savelabels{ii}];
    savefigure(fullfile(dataDir,'hierarchical_clustering',savename))
end

mask = ~cellfun(@(x) strcmp(x(1),'.'),genes);
Ydata = Ydata(mask,:); genes = genes(mask); foldChanges = foldChanges(mask);

%remove mitochondrial genes
geneI = find(cellfun(@(x) length(x) >= 3,genes));
geneslong = genes(geneI);
midxs = geneI(cellfun(@(x) strcmp(x(1:3),'MT-'),geneslong));
mask = true(size(Ydata,1),1); mask(midxs) = false;
Ydata = Ydata(mask,:); genes = genes(mask); foldChanges = foldChanges(mask);

%filter out low-expressed genes
meancpm = mean(Ydata,2);
if strcmp(dset,'time_cpm')
    thresh = -2.2;
elseif strcmp(dset,'time_cpm_lin')
    thresh = 2.5;
end

mask = meancpm > thresh;
Ydata = Ydata(mask,:);
genes = genes(mask);
foldChanges = foldChanges(mask);

%filter by fold change
thresh = std(foldChanges);

%save histogram of fold changes
figure('Position',figpos)
histogram(foldChanges)
cleanSubplot(fs); axis square; set(gca,'Box','off')
xline(-thresh,'LineWidth',2,'Color','k','LineStyle','--'); xline(thresh,'LineWidth',2,'Color','k','LineStyle','--');
xlabel('cumulative fold change'); ylabel('frequency')
xlim([-10,10])
savefigure(fullfile(dataDir,'hierarchical_clustering','foldChangeHistogram.png'))

mask = foldChanges < -thresh | foldChanges > thresh;
Ydata = Ydata(mask,:); genes = genes(mask); foldChanges = foldChanges(mask);
timefoldchange = foldChanges;

timegenes = genes; %keep for comparison with genes in the dose response
timemax = max(Ydata,[],2);

%preprocessing -> divide by the maximum
Ydata = Ydata./timemax;


%% build hierarchical tree
disp('getting distances')
disp(strcat("using ", dist, " distance"))
tic
Y = pdist(Ydata,dist);
toc
disp('doing linking')
disp(strcat("building clusters with ", method, " linkage"))
tic
Z = linkage(Y,method);
toc
disp('calculating cophenet')
tic
c = cophenet(Z,Y);

fprintf('cophenet = %.3g\n',c)
toc

%% split into determined number of clusters
k = 3;
colors = distinguishable_colors(k,'w');
cutoff = median([Z(end-k+1,3) Z(end-k+2,3)]);
Idx = cluster(Z,'Cutoff',cutoff,'criterion','distance');

figure('Position',figurePosition(560*1.3154,560))
ax1 = subplot_tight(1,2,1);
[~,~,outperm] = dendrogram(Z,0,'Orientation','left','ColorThreshold',cutoff);
invperm = 1:length(outperm);
invperm(outperm) = 1:length(outperm);

cleanSubplot(fs); axis square; xticklabels({}); yticklabels({})
ax1.YAxis.Visible = 'off'; ax1.XAxis.Visible = 'off';
ax1.PlotBoxAspectRatio = [0.25 1 1];
ax2 = subplot_tight(1,2,2);
imagesc(Ydata(outperm,:))
cleanSubplot; axis square; colormap turbo
yticklabels({}); xticklabels({}); ax2.Box = 'off';
ax2.YAxis.Visible = 'off'; ax2.XAxis.Visible = 'off';
set(ax2,'Ydir','normal')
set(ax1,'Position',[0.01,0.01,0.245,0.98])
set(ax2,'Position',[0.245,0.01,0.745,0.98])

sIdx = Idx(outperm);
diffs = zeros(size(sIdx)); diffs(2:end) = sIdx(2:end) - sIdx(1:end-1);
diffI = find(abs(diffs) > 0);
for ii = 1:length(diffI)
    yline(diffI(ii),'Color','w','LineWidth',3,'Alpha',1);
end

savename = strcat('treeandheatmap_',dset,'_',dist,'_',method);
savename = fullfile(dataDir,'hierarchical_clustering',savename);
savefigure(savename)

figure('Position',figurePosition(560*3,560))
subplot(1,3,1)
dendrogram(Z,'ColorThreshold',cutoff);
% [~, ~, outperm] = dendrogram(Z);
title('dendrogram')
cleanSubplot(fs); axis square
xticklabels({})

subplot(1,3,2)
histogram(Idx)
cleanSubplot(fs); axis square
set(gca,'XTick',1:k,'XTickLabel',cellstr(num2str((1:k)'))')
title('histogram')

subplot(1,3,3); hold on
for ki = 1:k
    if sum(Idx == ki) > 10
        plot(times,mean(Ydata(Idx == ki,:),1),'LineWidth',3,...
            'Color',colors(ki,:))
    end
end
cleanSubplot(fs); axis square
% xlabel('time (hours)'); ylabel('log_2cpm')
xlabel('time (hours)'); ylabel('cpm')
ylim([0,1])
% legend('Position',[0.915,0.6125,0.075,0.2])
lgd = legend; lgd.Position(1) = 0.915;
title('cluster centroids')

savename = strcat('hierarchical_',dset,'_',dist,'_',method,sprintf('_k%.2d.png',k));
savename = fullfile(dataDir,'hierarchical_clustering',savename);
savefigure(savename)

%% plot dendrogram and heatmap and show genes

nt = size(Ydata,2);
treeratio = 0.25;
iwidth = 0.98/(1 + 1/nt + treeratio);
width1 = iwidth*treeratio; width2 = iwidth*(1 + 1/nt);

handgenes = {'ISL1','SOX2','NANOG','TFAP2C','GATA3','HAND1','POU5F1'};

Is = cellfun(@(x) find(strcmp(genes,x)),handgenes,'UniformOutput',false);
mask = ~cellfun(@isempty,Is);
hgenes = handgenes(mask); Is = cell2mat(Is(mask));

ycoords = invperm(Is);

figure('Position',figurePosition(560/(iwidth/0.98),560))
ax1 = subplot_tight(1,3,1);
[~,~,outperm] = dendrogram(Z,0,'Orientation','left','ColorThreshold',cutoff);
invperm = 1:length(outperm);
invperm(outperm) = 1:length(outperm);

cleanSubplot(1.5*fs); axis square; xticklabels({}); yticklabels({})
ax1.YAxis.Visible = 'off'; ax1.XAxis.Visible = 'off';
ax1.PlotBoxAspectRatio = [0.25 1 1];
ax2 = subplot_tight(1,3,2);
imagesc(Ydata(outperm,:))
cleanSubplot(1.5*fs); axis square; colormap turbo
yticklabels({}); xticklabels({});  ax2.Box = 'off';
ax2.YAxis.Visible = 'off'; ax2.XAxis.Visible = 'off';
% xticklabels({'0','4','12','20','28','42'});
set(ax2,'Ydir','normal')

sIdx = Idx(outperm);
diffs = zeros(size(sIdx)); diffs(2:end) = sIdx(2:end) - sIdx(1:end-1);
diffI = find(abs(diffs) > 0);
for ii = 1:length(diffI)
    yline(diffI(ii),'Color','k','LineWidth',5,'Alpha',1,'LineStyle','-');
    yline(diffI(ii),'Color','w','LineWidth',3,'Alpha',1,'LineStyle','-');
end


ax3 = subplot_tight(1,3,3);

xlim([0,1]); ylim([1,size(Ydata,1)])
cleanSubplot(1.5*fs); yticklabels({}); xticklabels({}); ax3.Box = 'off';
ax3.YAxis.Visible = 'off'; ax3.XAxis.Visible = 'off';

set(ax1,'Position',[0.01,0.01,width1,0.98])
set(ax2,'Position',[width1+0.01,0.01,iwidth,0.98])
set(ax3,'Position',[width1+iwidth+0.01,0.01,iwidth/nt,0.98])

tplot = text(0.1*ones(size(ycoords)),ycoords,hgenes,'FontWeight','bold','FontSize',14);
hold on
plot([0;0.08],repmat(ycoords,2,1), 'LineWidth',2, 'Color', 'k');
for ii = 1:length(tplot), tplot(ii).Visible = 'off'; end

savename = strcat('tree_heatmap_genes_',dset,'_',dist,'_',method,sprintf('_k%.2d.png',k));
savename = fullfile(dataDir,'hierarchical_clustering',savename);
savefigure(savename)

for ii = 1:length(tplot), tplot(ii).Visible = 'on'; end

savename = [savename(1:end-4),'_withlabels.png'];
savefigure(savename)


%% SMAD4 levels (dataset with a single time point and multiple signal levels)
dset = 'dose_cpm_lin';
Ydata = M.(dset);

%remove genes with names beginning in the character '.' (novel transcripts
%and pseudogenes)
if strcmp(dset,'dose_cpm_lin') || strcmp(dset,'dose_cpm')
    genes = T.dose_cpm.Symbol;
elseif strcmp(dset,'dose_fpkm_lin') || strcmp(dset,'dose_fpkm')
    genes = T.dose_fpkm.Symbol;
end

%remove novel transcripts and pseudogenes
mask = ~cellfun(@(x) strcmp(x(1),'.'),genes);
Ydata = Ydata(mask,:); genes = genes(mask);
%remove mitochondrial genes
geneI = find(cellfun(@(x) length(x) >= 3,genes));
geneslong = genes(geneI);
midxs = geneI(cellfun(@(x) strcmp(x(1:3),'MT-'),geneslong));
mask = true(size(Ydata,1),1); mask(midxs) = false;
Ydata = Ydata(mask,:); genes = genes(mask);


Ydata = Ydata(:,1:end-1);

%only consider genes detected in both datasets
[C,ia,ib] = intersect(genes,timegenes);
Ydata = Ydata(ia,:)./timemax(ib);
genes = genes(ia);
meancpm = mean(Ydata,2);
ngene = size(Ydata,1);
tgenes = timegenes(ib); tfoldchange = timefoldchange(ib);
tcluster = Idx(ib);

if strcmp(dset(end-3:end),'_lin')
    foldChanges = log2(Ydata(:,1)./Ydata(:,5));
else
    foldChanges = Ydata(:,1) - Ydata(:,5);
end

% plot SOX2 vs SMAD4 integral with an example fit line
figure('Position',figurePosition(795,536)); hold on
x = SMAD4(1:end-1)'*4.5/42;
y = Ydata(strcmp(genes,'SOX2'),:)';
plot(x,y,'x','LineWidth',lw,'MarkerSize',15)

xinv = pinv([x,ones(size(x))]);
mb = xinv*y;
plot(x,mb(1)*x + mb(2),'LineWidth',lw,'LineStyle','--','Color',[0,0,0,0.65])

cleanSubplot(fs)%; axis square
xlabel('SMAD4 (N:C) integral'); ylabel('SOX2 (norm CPM)')
xlim([min(x),max(x)])

savename = ['correlationfitex_',dset,'_SOX2_v_int.png'];
savefigure(fullfile(dataDir,'dose_response',savename))

%% correlation with SMAD4
smad4 = SMAD4(1:end-1);

R = NaN(ngene,1);
for ii = 1:ngene
    r = corrcoef([Ydata(ii,:)',smad4']);
    R(ii) = r(1,2);
end

bthresh = 0.1; rthresh = 0.9;
%least squares fit of cpm to smad4 N:C level for each gene to calculate
%slopes scaled to maximum transcript levels in the time-series data
pA = pinv([smad4', ones(length(smad4),1)]);
mb = (pA*Ydata')';
b = mb(:,1);
bscale = b./max(Ydata,[],2);

%list of hand-picked genes to highlight on the scatterplot:
mygenes = {'ISL1','SOX2','NANOG','TFAP2C','GATA3','HAND1','POU5F1'};
showgenes = true;

TFmask = false(size(genes));
for ii = 1:length(genes)
    if sum(strcmp(genes{ii},TFs)) > 0
        TFmask(ii) = true;
    end
end


% figure('Position',figurePosition([795,894]))
figure('Position',figurePosition([795,536]))
scatter(R,b,10,repmat([0.5 0.5 0.5],length(genes),1),'filled'); hold on
scatter(R(TFmask),b(TFmask),15,repmat(lines(1),sum(TFmask),1),'filled')
cleanSubplot(fs)
xlabel('corr with SMAD4'); ylabel('slope (normalized)')
ylim([-0.75,0.75])
xline(-rthresh,'LineWidth',1.5); xline(rthresh,'LineWidth',1.5);
yline(bthresh,'LineWidth',1.5); yline(-bthresh,'LineWidth',1.5); %yline(0);
% set(gca,'Position',[0.1755 0.1745 0.7295 0.7505])

if showgenes && true
    hold on
    I = ismember(genes,mygenes);
    x = R(I); y = b(I);
    scatter(x,y,100,'r','filled')
    
    savename = strcat('corrvslope_',dset);
    savename = fullfile(dataDir,'dose_response',savename);
    savefigure(savename)
    
    text(x,y,genes(I))
end

savename = strcat('corrvslopeannotated_',dset);
savename = fullfile(dataDir,'dose_response',savename);
savefigure(savename)

%smaller figure with full axis limits
figure('Position',figpos)
scatter(R,b,'.')
cleanSubplot(fs); axis square
xlabel('corr with SMAD4 level'); ylabel('slope (normalized)')
% xline(-rthresh); xline(rthresh); yline(bthresh); yline(-bthresh);
ylim([min(b),max(b)])
savename = strcat('corrvslope_',dset,'_small');
savename = fullfile(dataDir,'dose_response',savename);
savefigure(savename)

%% make a table of genes of interest and export to a csv file

gene_name = tgenes;
cluster_id = tcluster;
slope = b;
correlation = R;

%label potential integrators
nidx = tcluster(strcmp(tgenes,'SOX2')); pidx = tcluster(strcmp(tgenes,'GATA3')); %<- indices of negatively and positively regulated clusters of interest
posgenes = R > rthresh & b > bthresh;
neggenes = R < -rthresh & b < -bthresh;

integrator = false(size(gene_name));
integrator((cluster_id == pidx) & posgenes) = true;
integrator((cluster_id == nidx) & neggenes) = true;

cluster_name = cell(size(gene_name));
for ii = 1:length(gene_name)
    if cluster_id(ii) == pidx
        cluster_name{ii} = 'early increasing';
    elseif cluster_id(ii) == nidx
        cluster_name{ii} = 'decreasing';
    else
        cluster_name{ii} = 'late increasing';
    end
end

CSVtable = table(gene_name, cluster_id, cluster_name, slope, correlation, integrator);

writename = 'integratorgenecandidates.csv';
writename = fullfile(dataDir,writename);
writetable(CSVtable,writename);




