clear; close all; clc

%% load signaling and gene expression data
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
dataDir = scriptPath;

dirs = fullfile(scriptPath,{'220609','220629','220812'});

ndirs = length(dirs);
signalData = cell(1,ndirs);
Sall = cell(1,ndirs);
X = cell(1,ndirs); Xnorm = cell(1,ndirs); Xs = cell(1,ndirs); Xnorms = cell(1,ndirs);
tvecs = cell(1,ndirs);
ms = cell(1,ndirs); ratios = cell(1,ndirs); pamnion = cell(1,ndirs); Rall = cell(1,ndirs);
mcon = cell(1,ndirs);
conditions = cell(1,ndirs);

allintegral = cell(1,ndirs); alllevel = cell(1,ndirs);
allduration = cell(1,ndirs); allL0 = cell(1,ndirs);

%load signaling and cell fate data
for di = 1:ndirs
    %signaling data
    sdata = load(fullfile(dirs{di},'allSignalData.mat'));
    signalData{di} = sdata.signalData;
    Sall{di} = sdata.Sall;
    tvecs{di} = sdata.tvec;
    %gene expression/cell fate data
    gdata = load(fullfile(dirs{di},'combinedExpressionData')); %should metadata also be included here?
    ms{di} = gdata.m;
    conditions{di} = gdata.conditions;
    
    gdata = load(fullfile(dirs{di},'ControlConditionsIF.mat'));
    mcon{di} = gdata.m;
end
ntime = min(cellfun(@(x) size(x,1),signalData));
tvec = tvecs{1}(1:ntime);

%normalize signaling data
maxconds = [1 1 9];
for di = 1:ndirs
    signalData{di} = signalData{di}(1:ntime,:);
    Sall{di} = Sall{di}(1:ntime,:,:);
    X{di} = cellfun(@(x) mean(x,'omitnan'),signalData{di});
    Xs{di} = cellfun(@(x) mean(x,'omitnan'),Sall{di});
    
    maxval = mean(X{di}(:,maxconds(di)));
    minval = mean(X{di}(tvec <= 0,:),'all');
    Xnorm{di} = (X{di} - minval)/(maxval - minval);
    Xnorms{di} = (Xs{di} - minval)/(maxval - minval);
end

%normalize expression ratio and get percent amnion
chans = [2 4; 2 4; 2 4]; %ISL1 and NANOG channels
tol = 0.001; %call bottom 0.1% of values baseline
threshs = NaN(1,ndirs); bws = NaN(1,ndirs); %threshold and bandwidth
for di = 1:ndirs
    M = cell2mat(mcon{di}(:));
    x = M(:,chans(di,1)); y = M(:,chans(di,2));
    [xl, yl] = axlims(x,y,tol,false);
    x = x - xl(1); x(x < 0) = 0; x = x + 1;
    y = y - yl(1); y(y < 0) = 0; y = y + 1;
    ratio = log10(x./y);
    rl = axlims(ratio,ratio,tol,false);
    % use k means to separate into ISL1+/-(equivalent to gmm/otsu thresholding in 1D)
    Idx = kmeans(ratio,2);
    t1 = min([max(ratio(Idx == 1)),max(ratio(Idx == 2))]);
    t2 = max([min(ratio(Idx == 1)),min(ratio(Idx == 2))]);
    thresh = 0.5*(t1 + t2);
    threshs(di) = thresh;
    %find peak locations for log(ISL1 / NANOG) density
    pts = linspace(rl(1),rl(2),100);
    [f,xi] = ksdensity(ratio);
    [~,locs,~,p] = findpeaks(f,xi); %get location and prominence of peaks
    %get the two peaks with greatest prominence (p)
    [~,I] = sort(p,'descend'); %order the peaks by prominence
    rpeaks = locs(I(1:2)); %get the two most prominent peaks 
    bw = 0.5*abs(rpeaks(1) - rpeaks(2));
    bws(di) = bw;
    rs = cell(size(ms{di}));
    ps = NaN(size(ms{di}));
    for ii = 1:size(ms{di},1)
        for jj = 1:size(ms{di},2)
            x = ms{di}{ii,jj}(:,chans(di,1));
            y = ms{di}{ii,jj}(:,chans(di,2));
            x = x - xl(1); x(x < 0) = 0; x = x + 1;
            y = y - yl(1); y(y < 0) = 0; y = y + 1;
            ratio = log10(x./y);
            ratio = (ratio - thresh)/bw;
            if numel(ratio) > 10
                rs{ii,jj} = ratio;
                ps(ii,jj) = 100*sum(ratio > 0)/numel(ratio);
            end
        end
    end
    ratios{di} = rs;
    pamnion{di} = ps;
    Rs = cell(size(rs,1),1);
    for ii = 1:size(rs,1)
        Rs{ii} = cell2mat(rs(ii,:)');
    end
    Rall{di} = Rs;
end

%% options for all plots
opts = struct(...
    'fs',       28,... %font size
    'lfs',      20,... %legend font size
    'figpos',   figurePosition([560 560]),... %figure size
    'lw',       3,... %line width
    'yscale',   0.95,... %scaling factor for the axes in y
    'xscale',   1,...
    'savedir',  fullfile(dataDir,'figures')); %scaling factor for the axes in x

defaultcolors = [0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125;...
    0.494 0.184 0.556; 0.466 0.674 0.188];

%sigmoid model used for finding thresholds
ft = fittype('100/(1 + exp(-lambda*(t-xmid)))','indep','t');
% coeffnames(ft)
startPoint = [25,0.5];

%% plot signaling data - 220629 dataset
opts.colororder = 'reverse';

di = 2;
xdata = Xnorm{di};

conds = {1:8,9:16};
opts.titles = {{'signaling:','full duration (42 hr)'},{'signaling:','reduced duration (32 hr)'}};
opts.savenames = {'signaling_fullduration','signaling_reducedduration'};
plotSignaling(xdata,tvec,conds,opts)

%% plot expression ratio data - 220629 dataset
opts.titles = {{'cell fate:','full duration'},{'cell fate:','reduced duration'}};
legstr = conditions{di}(conds{1})';
legstr = strrep(legstr,'LDN',''); legstr = strrep(legstr,' 42hr','');
opts.legstr = legstr;
% opts.legtitle = 'BMPRi (nM)';
opts.legtitle = {'BMPRi','(nM)'};
opts.savenames = {'fate_fullduration','fate_reducedduration'};

plotRatioDistribution(Rall{di},conds,opts)

opts.legstr = '';
opts.savenames = {'fate_fullduration_nolegend','fate_reducedduration_nolegend'};
plotRatioDistribution(Rall{di},conds,opts)

%% level thresholds - 220629 dataset
%get levels before and after LDN from the data

xdata = Xnorms{di};
levels = NaN(size(xdata,2),size(xdata,3));
L0s = NaN(size(xdata,2),size(xdata,3));
integrals = NaN(size(xdata,2),size(xdata,3));
durations = NaN(size(xdata,2),size(xdata,3));
tmaxes = NaN(size(xdata,2),1);
tmaxes(conds{1}) = 42; tmaxes(conds{2}) = 32;

duridxs = conds{2}(1:end-1);
for ii = 1:size(xdata,2)
    tmask = tvec > 0 & tvec <= tmaxes(ii); %time points over which to do the calculation
    for jj = 1:size(xdata,3)
        levels(ii,jj) = mean(xdata(tmask,ii,jj),'omitnan');
        integrals(ii,jj) = mean(xdata(tvec > 0,ii,jj),'omitnan');
        if ismember(ii,duridxs)
            durations(ii,jj) = getduration(tvec(tvec>0),xdata(tvec>0,ii,jj));
        else
            durations(ii,jj) = 42;
        end
        L0s(ii,jj) = mean(xdata(tvec > durations(ii,jj)+1,ii,jj));
    end
end

allintegral{di} = integrals; allduration{di} = durations;
alllevel{di} = levels; allL0{di} = L0s;

%level thresholds
levelthresh = NaN(length(conds),1);
lthresh = NaN(length(conds),1);% lambdas = NaN(length(conds),1);
lerr = NaN(length(conds),1);
figure('Position',opts.figpos); hold on
p = gobjects(length(conds),1);
for ii = 1:length(conds)
    xall = levels(conds{ii},:);
    yall = pamnion{di}(conds{ii},:);
    
    xdata = mean(xall,2,'omitnan');
    xerr = std(xall,0,2,'omitnan');
    ydata = mean(yall,2,'omitnan');
    yerr = std(yall,0,2,'omitnan');
    
    errorbar(xdata,ydata,yerr,yerr,xerr,xerr,'o','LineWidth',opts.lw,...
        'Color',defaultcolors(ii,:))
    levelthresh(ii) = interp1(ydata,xdata,50);
    
    %threshold by interpolation
    levelthresh(ii) = interp1(ydata,xdata,50);
    %threshold by sigmoid fit
    xall = xall(:); yall = yall(:);
    nanmask = any(isnan([xall,yall]),2);
    xall = xall(~nanmask); yall = yall(~nanmask);
    mdl = fit(xall,yall,ft,'StartPoint',startPoint);
    lthresh(ii) = mdl.xmid; lambda = mdl.lambda;
    test = confint(mdl);
    lerr(ii) = test(2,2) - lthresh(ii);
    xx = linspace(min(xdata),max(xdata),100);
    p(ii) = plot(xx, 100./(1 + exp(-lambda*(xx - lthresh(ii)))),'LineWidth',opts.lw,...
        'Color',[defaultcolors(ii,:),1]);
end

for ii = 1:length(conds)
    xline(lthresh(ii),'LineWidth',2,'LineStyle','-','Color',[defaultcolors(ii,:),0.5]);
end
cleanSubplot(opts.fs); axis square
legend(p,'full duration','reduced duration','Location','northwest','FontSize',opts.lfs)
xlabel('SMAD4 level (au)'); ylabel('percent amnion')
title({'','level threshold'})
ylim([0,100]); xlim([0,1.1])

scaleaxes(opts.xscale,opts.yscale)
savefigure(fullfile(opts.savedir,'level_thresholds_220629'))

colidxs = [1 2];%colidxs = [3 4];
%integral thresholds
intthresh = NaN(length(conds),1);
ithresh = NaN(length(conds),1);
ierr = NaN(length(conds),1);
figure('Position',opts.figpos); hold on
p = gobjects(length(conds),1);
for ii = 1:length(conds)
    xall = integrals(conds{ii},:);
    yall = pamnion{di}(conds{ii},:);
    
    xdata = mean(xall,2,'omitnan');
    xerr = std(xall,0,2,'omitnan');
    ydata = mean(yall,2,'omitnan');
    yerr = std(yall,0,2,'omitnan');
    errorbar(xdata,ydata,yerr,yerr,xerr,xerr,'o','LineWidth',opts.lw,...
        'Color',defaultcolors(colidxs(ii),:))
    intthresh(ii) = interp1(ydata,xdata,50);
    %threshold by sigmoid fit
    xall = xall(:); yall = yall(:);
    nanmask = any(isnan([xall,yall]),2);
    xall = xall(~nanmask); yall = yall(~nanmask);
    mdl = fit(xall,yall,ft,'StartPoint',startPoint);
    ithresh(ii) = mdl.xmid; lambda = mdl.lambda;
    test = confint(mdl);
    ierr(ii) = test(2,2) - ithresh(ii);
    xx = linspace(min(xdata),max(xdata),100);
    p(ii) = plot(xx, 100./(1 + exp(-lambda*(xx - ithresh(ii)))),'LineWidth',opts.lw,...
        'Color',[defaultcolors(colidxs(ii),:),1]);
end

for ii = 1:length(conds)
    xline(ithresh(ii),'LineWidth',2,'LineStyle','-',...
        'Color',[defaultcolors(colidxs(ii),:),0.5]);
end
cleanSubplot(opts.fs); axis square
xlabel('SMAD4 integral (au)'); ylabel('percent amnion')
title({'','integral threshold'})
ylim([0,100]); xlim([0,1.1]);

scaleaxes(opts.xscale,opts.yscale)
savefigure(fullfile(opts.savedir,'integral_thresholds_220629'))

fprintf('level threshold 1: %.3g +/- %.2g\n',lthresh(1),lerr(1))
fprintf('level threshold 2: %.3g +/- %.2g\n',lthresh(2),lerr(2))
fprintf('integral threshold 1: %.3g +/- %.2g\n',ithresh(1),ierr(1))
fprintf('integral threshold 2: %.3g +/- %.2g\n',ithresh(2),ierr(2))

%% plot signaling data - 220812 dataset
opts.colororder = 'normal';

di = 3;
xdata = Xnorm{di};
xdata(:,18) = Xnorm{1}(:,4);

opts.savenames = {'signaling_highlevel','signaling_reducedlevel'};

conds = {1:9,10:18};
opts.titles = {{'signaling:','high level (0 nM)'},{'signaling:','reduced level (10 nM)'}};
plotSignaling(xdata,tvec,conds,opts)

%% plot expression ratio data - 220812 dataset
opts.titles = {{'cell fate:','high level'},{'cell fate:','reduced level'}};
ratio = Rall{di};
ratio{18} = Rall{1}{4};

legstr = conditions{di}(conds{1})';
legstr = strrep(legstr,'LDN0,',''); legstr = strrep(legstr,'hr','');
opts.legstr = legstr;
opts.legtitle = 'hours';

opts.savenames = {'fate_highlevel','fate_reducedlevel'};
plotRatioDistribution(ratio,conds,opts)

opts.legstr = '';
opts.savenames = {'fate_highlevel_nolegend','fate_reducedlevel_nolegend'};
plotRatioDistribution(ratio,conds,opts)


%% duration thresholds - 220812 dataset
%for durations other than 0 and 42hr, get signaling duration from the data

xdata = Xnorms{di};
xdata(:,18,:) = repmat(Xnorms{1}(:,4,:),1,1,2);

fatedata = pamnion{di};
fatedata(18,:) = repmat(pamnion{1}(4,2:3),1,4);

durations = NaN(size(xdata,2),size(xdata,3));
L0s = NaN(size(xdata,2),size(xdata,3));
integrals = NaN(size(xdata,2),size(xdata,3));
levels = NaN(size(xdata,2),size(xdata,3));

%only calculate duration for intermediates, not conditions with 0 or 42hr
%duration
times = tvec(tvec > 0);
for ii = 1:size(xdata,2)
    for jj = 1:size(xdata,3)
        ydata = xdata(tvec > 0,ii,jj);
        if sum(isnan(ydata)) < 10
            if ismember(ii,[1,10]) %0 hour conditions
                durations(ii,jj) = 0;
                tmax = 42;
            elseif ismember(ii,[9,18]) %42 hour conditions
                durations(ii,jj) = 42;
                tmax = 42;
            else
                durations(ii,jj) = getduration(times,ydata);
                tmax = durations(ii,jj);
            end
            integrals(ii,jj) = mean(xdata(tvec > 0,ii,jj),'omitnan');
            levels(ii,jj) = mean(xdata(tvec > 0 & tvec < tmax - 0.5,ii,jj),'omitnan');
            L0s(ii,jj) = mean(xdata(tvec > tmax+2,ii,jj));
        end
    end
end

allintegral{di} = integrals; allduration{di} = durations;
alllevel{di} = levels; allL0{di} = L0s;

%duration thresholds
durthresh = NaN(length(conds),1); dthresh = NaN(length(conds),1);
derr = NaN(length(conds),1);
figure('Position',opts.figpos); hold on
p = gobjects(length(conds),1);
startPoint = [1 30];
for ii = 1:length(conds)
    xall = durations(conds{ii},:);
    yall = fatedata(conds{ii},:);
    
    xdata = mean(xall,2,'omitnan');
    xerr = std(xall,0,2,'omitnan');
    ydata = mean(yall,2,'omitnan');
    yerr = std(yall,0,2,'omitnan');
    errorbar(xdata,ydata,yerr,yerr,xerr,xerr,'o','LineWidth',opts.lw,...
        'Color',defaultcolors(ii,:))
    %threshhold by interpolation
    durthresh(ii) = interp1(ydata,xdata,50);
    
    %threshold by sigmoid fit
    xall = xdata; yall = ydata;
    mdl = fit(xall,yall,ft,'StartPoint',startPoint);
    dthresh(ii) = mdl.xmid; lambda = mdl.lambda;
    test = confint(mdl);
    derr(ii) = test(2,2) - dthresh(ii);
    xx = linspace(min(xdata),max(xdata),100);
    p(ii) = plot(xx, 100./(1 + exp(-lambda*(xx - dthresh(ii)))),'LineWidth',opts.lw,...
        'Color',[defaultcolors(ii,:),1]);
end

for ii = 1:length(conds)
    xline(dthresh(ii),'LineWidth',2,'LineStyle','-','Color',[defaultcolors(ii,:),0.5]);
end
cleanSubplot(opts.fs); axis square
legend(p,'high level','reduced level','Location','northwest','FontSize',opts.lfs)
xlabel('SMAD4 duration (hr)'); ylabel('percent amnion')
title({'','duration threshold'})
ylim([0,100]); xlim([0,42.5]);

scaleaxes(opts.xscale,opts.yscale)
savefigure(fullfile(opts.savedir,'duration_thresholds_220812'))

%integral thresholds
colidxs = [1 2];
intthresh = NaN(length(conds),1); ithresh = NaN(length(conds),1);
ierr = NaN(length(conds),1);
figure('Position',opts.figpos); hold on
startPoint = [25 0.7];
for ii = 1:length(conds)
    xall = integrals(conds{ii},:);
    yall = fatedata(conds{ii},:);
    
    xdata = mean(xall,2,'omitnan');
    xerr = std(xall,0,2,'omitnan');
    ydata = mean(yall,2,'omitnan');
    yerr = std(yall,0,2,'omitnan');
    errorbar(xdata,ydata,yerr,yerr,xerr,xerr,'o','LineWidth',opts.lw,...
        'Color',defaultcolors(colidxs(ii),:))
    %threshhold by interpolation
    intthresh(ii) = interp1(ydata,xdata,50);
    
    %threshold by sigmoid fit
    xall = xdata; yall = ydata;
    mdl = fit(xall,yall,ft,'StartPoint',startPoint);
    ithresh(ii) = mdl.xmid; lambda = mdl.lambda;
    test = confint(mdl);
    ierr(ii) = test(2,2) - ithresh(ii);
    xx = linspace(min(xdata),max(xdata),100);
    p(ii) = plot(xx, 100./(1 + exp(-lambda*(xx - ithresh(ii)))),'LineWidth',opts.lw,...
        'Color',[defaultcolors(colidxs(ii),:),1]);
end

for ii = 1:length(conds)
    xline(ithresh(ii),'LineWidth',2,'LineStyle','-',...
        'Color',[defaultcolors(colidxs(ii),:),0.5]);
end
cleanSubplot(opts.fs); axis square
xlabel('SMAD4 integral (au)'); ylabel('percent amnion')
title({'','integral threshold'})
ylim([0,100]); xlim([0,1.1]);

scaleaxes(opts.xscale,opts.yscale)
savefigure(fullfile(opts.savedir,'integral_thresholds_220812'))

fprintf('duration threshold 1: %.3g +/- %.2g\n',dthresh(1),derr(1))
fprintf('duration threshold 2: %.3g +/- %.2g\n',dthresh(2),derr(2))
fprintf('integral threshold 1: %.3g +/- %.2g\n',ithresh(1),ierr(1))
fprintf('integral threshold 2: %.3g +/- %.2g\n',ithresh(2),ierr(2))

%% local functions

function plotSignaling(xdata,tvec,conds,opts)

if ~isfield(opts,'xscale')
    xscale = 1;
else
    xscale = opts.xscale; 
end

if ~isfield(opts,'yscale')
    yscale = 1;
else
    yscale = opts.yscale; 
end

if ~isfield(opts,'savedir')
    savedir = [];
else
    savedir = opts.savedir; 
end

for ii = 1:length(conds)
    colors = turbo(length(conds{ii}));
    if strcmp(opts.colororder,'reverse')
        colors = colors(end:-1:1,:);
    end
    
    figure('Position',opts.figpos)
    hold on
    for cidx = 1:length(conds{ii})
        v = xdata(:,conds{ii}(cidx));
        plot(tvec,v,'LineWidth',opts.lw,'Color',colors(cidx,:))
    end
    hold off
    cleanSubplot(opts.fs); axis square
    xlim([min(tvec)-0.5,max(tvec)+0.5])
    ylim([-0.05,1.3])
    ylabel('SMAD4 (N:C)')
    xlabel('time (hr)')%xlabel('time (hours)')
    title(opts.titles{ii})
    
    %adjust scale so everything fits in the plot
    scaleaxes(xscale,yscale)
    savefigure(fullfile(savedir,opts.savenames{ii}))
end

end


function plotRatioDistribution(ratio,conds,opts)

if ~isfield(opts,'xscale')
    xscale = 1;
else
    xscale = opts.xscale; 
end

if ~isfield(opts,'yscale')
    yscale = 1;
else
    yscale = opts.yscale; 
end

if ~isfield(opts,'savedir')
    savedir = [];
else
    savedir = opts.savedir; 
end

for ii = 1:length(conds)
    colors = turbo(length(conds{ii}));
    if strcmp(opts.colororder,'reverse')
        colors = colors(end:-1:1,:);
    end
    
    figure('Position',opts.figpos)
    hold on
    for cidx = 1:length(conds{ii})
        v = ratio{conds{ii}(cidx)};
        [f,xi] = ksdensity(v);
        plot(xi,f,'LineWidth',opts.lw,'Color',colors(cidx,:))
    end
    hold off
    cleanSubplot(opts.fs); axis square
    xlim([-2,2])
    ylabel('density estimate')
    xlabel('log(ISL1 / NANOG)') 
    title(opts.titles{ii})
    
    xline(0,'LineStyle','--','LineWidth',2);
    
    if isfield(opts,'legstr') && ~isempty(opts.legstr)
        lgd = legend(opts.legstr,'FontSize',opts.lfs);
        if isfield(opts,'legtitle') && ~isempty(opts.legtitle)
            title(lgd,opts.legtitle,'FontSize',opts.lfs)
        end
    end
    
    %adjust scale so everything fits in the plot
    scaleaxes(xscale,yscale)
    %add legend, adjust location
    if isfield(opts,'legstr') && ~isempty(opts.legstr)
        disp('paused: adjust legend position')
        pause
        fprintf('\n')
    end
    savefigure(fullfile(savedir,opts.savenames{ii}))
end


end


function scaleaxes(xscale,yscale)

pos = get(gca,'Position');
pos(2) = pos(2) + (1 - yscale)*pos(4);
pos(4) = yscale*pos(4);
pos(1) = pos(1) - (1 - xscale)*pos(3);
pos(3) = xscale*pos(3);
set(gca,'Position',pos)

end

function duration = getduration(times,ydata)

diffs = NaN(size(times));
for ti = 2:length(times)-1
    diffs(ti) = mean(ydata(times<times(ti)),'omitnan') -...
        mean(ydata(times>times(ti)),'omitnan');
end
[~,I] = max(diffs,[],'omitnan');
duration = times(I);

end





