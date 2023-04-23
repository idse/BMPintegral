clear; close all; clc

restoredefaultpath
rehash toolboxcache

scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(genpath(scriptPath));
addpath(genpath('/Users/idse/repos/HeemskerkLabCode'))

%% setup

%graphical options
fs = 24; %font size
wh = [560 560]; figpos = figurePosition(wh); %figure width and height -> figure position
lw = 3; %line width
tol = 0.01; %tolerance

scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
dataDir = fullfile(scriptPath,'matfiles');
savedir = fullfile(dataDir,'figures');

%dsetlabel = {'210801_histories','210827_histories','210827_corrected', '220928_histories'};%,'210827_corrected_only_updated','211115_histories','220111_histories'};
dsetlabel = {'210827_histories_nomedfilt','210827_histories_corrected_nomedfilt','220928_histories_nomedfilt','210827_histories_combined_nomedfilt'};

hist_raw = {};
hist_magic = {};
hist_dmagic = {};
hist_dmagic_hist = {};
hist_dmagic_fate = {};
hist_scramble_dmagic = {};
tvec = {};
fm = {};
fm_dmagic = {};
fm_scramble_dmagic = {};
fm_dmagic_hist = {};
fm_dmagic_fate = {};
p = {};
conditions = {};
channelLabels = {};

for i = 1:numel(dsetlabel) 

    s = load(fullfile(dataDir,[dsetlabel{i} '.mat']));

    field = 'NCratio'; channel = 1;
    X = s.Mats.(field)(:,:,channel+1);
    hist_raw{i} = X;

    fm{i} = s.Mats.fateMarkers;
    % remove DAPI for consistency with double magic, 
    % but not necessary for normalized markers where it is already removed
    fm{i} = fm{i}(2:end,:); 
    channelLabels{i} = s.fixedMeta.channelLabel(2:end);
    
%     fm{i} = s.normalizedMarkers;
%     channelLabels{i} = s.normalizedChannels;
    
    ntime = s.liveMeta.nTime;
    timeres = strsplit(s.liveMeta.timeInterval,' ');
    timeres = str2double(timeres{1});
    tscale = timeres/60;
    tvec{i} = (0:ntime-1)'*tscale; %shift to treatment time

    %position (condition) from which each cell came signaling data
    p{i} = s.Mats.positionIdx;

    conditions{i} = s.fixedMeta.conditions;
    
    hist_magic{i} = readmatrix([dsetlabel{i} '_magic.csv']);
    hist_magic{i} = hist_magic{i}(2:end,2:end)';
    
    hist_dmagic{i} = readmatrix([dsetlabel{i} '_dmagic_separate.csv']);
    hist_dmagic{i} = hist_dmagic{i}(2:end,2:end)';
    fm_dmagic{i} = readmatrix([dsetlabel{i} '_fm_dmagic_separate.csv']);
    fm_dmagic{i} = fm_dmagic{i}(2:end,2:end)';
    
    hist_dmagic_scramble{i} = readmatrix([dsetlabel{i} '_scramble_dmagic.csv']);
    hist_dmagic_scramble{i} = hist_dmagic_scramble{i}(2:end,2:end)';
    fm_dmagic_scramble{i} = readmatrix([dsetlabel{i} '_fm_scramble_dmagic.csv']);
    fm_dmagic_scramble{i} = fm_dmagic_scramble{i}(2:end,2:end)';
    
    hist_dmagic_fate{i} = readmatrix([dsetlabel{i} '_dmagic_fate.csv']);
    hist_dmagic_fate{i} = hist_dmagic_fate{i}(2:end,2:end)';
    fm_dmagic_fate{i} = readmatrix([dsetlabel{i} '_fm_dmagic_fate.csv']);
    fm_dmagic_fate{i} = fm_dmagic_fate{i}(2:end,2:end)';
    
    hist_dmagic_hist{i} = readmatrix([dsetlabel{i} '_dmagic_hist.csv']);
    hist_dmagic_hist{i} = hist_dmagic_hist{i}(2:end,2:end)';
    fm_dmagic_hist{i} = readmatrix([dsetlabel{i} '_fm_dmagic_hist.csv']);
    fm_dmagic_hist{i} = fm_dmagic_hist{i}(2:end,2:end)';
end


% redefine 220928 conditions which mean position 
di = 3;
p{di}(p{di} == 2) = 1;
p{di}(p{di} == 3) = 1;
p{di}(p{di} == 4) = 2;
p{di}(p{di} == 5) = 2;
p{di}(p{di} == 6) = 2;
conditions{di} = {[conditions{di}{1}  '_well1'], [conditions{di}{1}  '_well2']};

% redefine 220928 channel labels (ONLY FOR NOT NORMALIZED FATE)
channelLabels{di}{4} = 'DAPI2';
channelLabels{di}{8} = 'DAPI3';
channelLabels{di}{11} = 'DAPI4';
channelLabels{di}{12} = 'TFAP2C-2';
channelLabels{di}{13} = 'GATA3-2';
channelLabels{di}


%% signaling density plot over time from all cells (not tracked)

highSig = {};
lowSig = {};

for di = 3
    
    dsetlabel{di}
    lt = load([dsetlabel{di}(1:6) '_lt.mat']);
    lt = lt.lt;
    positions = lt.live_position;

    treatmentTime = 2;
    s = strsplit(lt.liveMeta.timeInterval,'min');
    ntime = positions(1).nTime;
    %tvec = ((0:ntime-1) - treatmentTime)*str2double(s{1})/60;
    npos = length(positions);

    NCRmin = 0;
    NCRmax = 1.5;
    meanSig = NaN(1,ntime); medSig = NaN(1,ntime);
    meanSigFromf = NaN(1,ntime);
    for ti = 1:ntime
        %ti
        data = cell(npos,1);
        for ii = 1:npos
            data{ii} = positions(ii).cellData(ti).NCratio(:,2);
            goodones = positions(ii).cellData(ti).labels == 1;
            data{ii} = data{ii}(goodones);
            %sum(isnan(positions(ii).cellData(ti).NCratio(:,2)))
        end
        data = cell2mat(data);
        data = data(data>NCRmin & data<NCRmax);
        meanSig(ti) = nanmean(data);
        medSig(ti) = nanmedian(data);
    end
    % normalize to mean of the tracked histories (vs mean of all cells)
    highSig{di} = mean(mean(hist_raw{di}(6:66,:)));% mean(meanSig(6:66));
    lowSig{di} =  mean(mean(hist_raw{di}(1:2,:)));%mean(meanSig(1:4));

    npts = 300;
    NCRmin = (NCRmin - lowSig{di})/(highSig{di}-lowSig{di});%.3;
    NCRmax = (NCRmax - lowSig{di})/(highSig{di}-lowSig{di});
    pts = linspace(NCRmin,NCRmax,npts);

    meanHist = mean((hist_raw{di}-lowSig{di})/(highSig{di}-lowSig{di}),2);

    X = NaN(npts,ntime); T = repmat(tvec{di}',npts,1); S = repmat(pts',1,ntime);
    for ti = 1:ntime
        %ti
        data = cell(npos,1);
        for ii = 1:npos
            data{ii} = positions(ii).cellData(ti).NCratio(:,2);
            goodones = positions(ii).cellData(ti).labels == 1;
            data{ii} = data{ii}(goodones);
            %sum(isnan(positions(ii).cellData(ti).NCratio(:,2)))
        end
        data = cell2mat(data);
        data = data(data>NCRmin & data<NCRmax);
        data = (data - lowSig{di})/(highSig{di}-lowSig{di});
        [f,~] = ksdensity(data,pts); %,'Support','positive'
        X(:,ti) = f;
        meanSig(ti) = nanmean(data);
        medSig(ti) = nanmedian(data);
        meanSigFromf(ti) = sum(pts.*f./sum(f));
    end

    figure
    surf(T,S,X,'LineStyle','none')
    fs = 26;
    cleanSubplot(fs); view(2); axis square; colormap turbo
    xlim(tvec{i}([1,end])); ylim(pts([1,end]))
    xlabel('time (hr)'); ylabel('SMAD4 (N:C)')
    %title(dsetlabel{i}(1:6))
    ylim([-0.5 2.5]);
    title('signaling distribution')
    hold on
    plot3(tvec{di}, meanHist, tvec{di}*0 + max(X(:)),'LineWidth',3,'Color','k')
    legend({'probability', 'mean'})
    plot3(tvec{i}, meanSig, tvec{i}*0 + max(X(:)),'LineWidth',3,'Color','k')
    %plot3(tvec{i}, medSig, tvec{i}*0 + max(X(:)),'LineWidth',3,'Color','r')
    %plot3(tvec{i}, meanSigFromf, tvec{i}*0 + max(X(:)),'LineWidth',3,'Color','b')
    hold off
    savefigure(fullfile(savedir,[dsetlabel{di}(1:6) '_distributionVtime']));

end

%% signaling density plot over time from histories

highSig = {};
lowSig = {};

for di = 1:4
    
    highSig{di} = mean(mean(hist_raw{di}(6:66,:)));% mean(meanSig(6:66));
    lowSig{di} =  mean(mean(hist_raw{di}(1:2,:)));%mean(meanSig(1:4));

	npts = 100;
    pts = linspace(-0.5,2.5,npts);
    ntime = size(hist_raw{di},1);
    
    X = NaN(npts,ntime); 
    T = repmat(tvec{di}',npts,1); 
    S = repmat(pts',1,ntime);
    
    meanSig = NaN(1,ntime);
    
    for ti = 1:ntime
        
        data = hist_raw{di}(ti,:);
        data = (data - lowSig{di})/(highSig{di}-lowSig{di});
        [f,~] = ksdensity(data,pts); %,'Support','positive'
        X(:,ti) = f;
        
        meanSig(ti) = nanmean(data);
        medSig(ti) = nanmedian(data);
        meanSigFromf(ti) = sum(pts.*f./sum(f));
    end

    figure
    surf(T,S,X,'LineStyle','none')
    fs = 26;
    cleanSubplot(fs); view(2); axis square; 
    colormap(flipud(brewermap([],"RdBu")));
    caxis([0 1.7])
    %colormap turbo
    xlim(tvec{di}([1,end])); ylim(pts([1,end]))
    xlabel('time (hr)'); ylabel('SMAD4 (N:C)')
    %title(dsetlabel{i}(1:6))
    ylim([-0.5 2.5]);
    title(' ')
    hold on
    plot3(tvec{di}, meanSig, tvec{di}*0 + max(X(:)),'LineWidth',3,'Color','k')
    hold off
    grid off
    box on
    set(gca,'layer','top');
    legend({'probability', 'mean'})
    savefigure(fullfile(savedir,[dsetlabel{di} '_distributionVtimeFromHist']));
    close;
end

%% make fate marker correlation heatmap

%channels in order for heatmap

di = 3;
dsetlabel{di}
lt = load([dsetlabel{di}(1:6) '_lt.mat']);
lt = lt.lt;

liveMeta = lt.liveMeta;
fixedMeta = lt.fixedMeta;

channelLabel = fixedMeta.channelLabel;
npos = length(lt.live_position);

chans = [4 3 6 13 14 2 8];
nc = length(chans);
cls = channelLabel(chans);

m = cell(npos,1);
for ii = 1:npos
    %get expression of non-dividing SMAD4 cells
    mask = lt.fixed_position(ii).cellData.labels == 1;
    mask = mask & lt.mapped_idxs{ii}' > 0;
    m{ii} = lt.fixed_position(ii).cellData.nucLevel(mask,chans);
end

M = cell2mat(m);

dtype = 'logt';

if strcmp(dtype,'raw')
    %raw data
    X = M;
elseif strcmp(dtype,'logt')
    %log transformed data
    M = M - min(M,[],1) + 1;
    M = M./mean(M,1)*100; %scale so offset by one doesn't distort the distribution
    X = log10(M+1);
end

%z normalize
X = (X - mean(X,1))./std(X,0,1);

idxs = 1:size(X,2);
figure('Position',figurePosition([600 560]))
imagesc(corrcoef(X(:,idxs)),[-1 1])
cleanSubplot(30); colormap(flipud(brewermap([],"RdBu")))
set(gca, 'YTick', 1:length(idxs), 'YTickLabel', cls(idxs));
set(gca, 'XTick', 1:length(idxs), 'XTickLabel', cls(idxs));
xtickangle(gca,45)
colorbar('Ticks',[-1 0 1]);
title(' ')
%saveas(gcf,fullfile(savedir,'pairwise_correlations.png'))
close;
    
%% fit a 2D GMM for visualization

k = 2; %number of components
suffix = sprintf('_k%d',k);

%channels for which to visualize the gaussian fit (6 = isl1, 1 = nanog)
chan = [6 1];
X = log10(M+1);
xyz = X(:, chan);
cl = channelLabels{di}(chan);

options = statset('Display','final');
gm = fitgmdist(xyz,k,'Options',options);

figure
scatter(xyz(:,1),xyz(:,2),10,'k','filled')
hold on
gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(gm,[x0,y0]),x,y);
%XL = [-2 2]; YL = [-2 2];
[XL,YL] = axlims(xyz(:,1),xyz(:,2),0.01,false);
XL = [0.9 2.7];
YL = [0.4 2.9;];
h = fcontour(gmPDF,[XL, YL],'LineWidth',2);
fs = 26;
cleanSubplot(fs); axis square
xlabel(['log( ' cls{chan(1)} ' )']); ylabel(['log( ' cls{chan(2)} ' )'])
%mu = gm.mu;
%scatter(mu(:,1),mu(:,2),50,'r','filled')
% legend('data','gmm','centroids','Location','southwest')
xlim(XL); ylim(YL)
saveas(gcf,fullfile(savedir,['2d_gmm_fit',suffix,'.png']))


%% plot some random (or not random) traces 

di = 3;
figure,
for filtertype = 1:2
    
    dsetlabel{di}
    nhist = 50;
    %ridx = 1 + round(rand(nhist,1)*(size(hist{i},2))-1);
    ridx = 1:nhist;
    if filtertype==1
        histsample = hist_raw{di}(:,ridx);
        suffix = '';
        titlestr = ' ';%'individual histories';
    elseif filtertype==2
        histsample = hist_magic{di}(:,ridx);
        suffix = '_magic';
        titlestr = ' ';%'denoised histories';
    end
    histsample = (histsample - lowSig{di})/(highSig{di}-lowSig{di});
    colors = lines(nhist);
    clf
    hold on
    for ii = 1:nhist
        plot(tvec{di}, histsample(:,ii), 'LineWidth',1, 'Color', colors(ii,:))
    end
    plot(tvec{di}, meanSig, 'LineWidth',3,'Color','k')
    hold off
    xlim([tvec{i}(1), tvec{di}(end)])
    ylim([-0.5 2.5]);
    text(2,-0.25,['N=' num2str(nhist) '/' num2str(size(hist_raw{di},2))],'FontSize',fs)
    %ylim([0.4 1.4]);
    fs = 26;
    lw = 2;
    cleanSubplot(fs,lw);  axis square;
    title(titlestr)
    xlabel('time (hr)'); 
    yl = ylabel('SMAD4 (N:C)');
    set(yl, 'Units', 'Normalized', 'Position', [-0.2, 0.5, 0]);
    savefigure(fullfile(savedir,[dsetlabel{di}(1:6) '_' num2str(nhist) 'randomhistories' suffix]));
end


%% show fit

X = hist_raw{di};
X = (X - lowSig{di})/(highSig{di}-lowSig{di});

X_mag = hist_magic{di};
X_mag = (X_mag - lowSig{di})/(highSig{di}-lowSig{di});

visualizefits = 'on';
nhists = size(X,2);

Tm = zeros(nhists,1); Us = Tm; Ls = Tm; Ds = Tm;
            
startPoint = [0.5 1 26 1]; %starting values from which to optimize the fit
colors = lines(2);
tic
magictoo = true;
figure
for ii = 10%89%1:nhists % 89, 10
    clf
    v = X(:,ii);
    tic
    t = tvec{di}(4:end); % exclude time before BMP
    % x(1) = high; x(2) = high - low; x(3) = duration; x(4) = sharpness of
    % drop
    fun = @(x) x(1) - x(2)./(1 + exp(-(t-x(3))/x(4)));
    fundiff = @(x) fun(x) - X(4:end,ii);
    fundiff_mag = @(x) fun(x) - X_mag(4:end,ii);
    x0 = startPoint;
    lb = [0 0 0 0];
    ub = [2 2 42 42];
    fitparam = lsqnonlin(fundiff, x0, lb, ub);
    fitparam_mag = lsqnonlin(fundiff_mag, x0, lb, ub);
    Us(ii) = fitparam(1); Ls(ii) = fitparam(1)- fitparam(2); Tm(ii) = fitparam(3); 
    
    toc
    if strcmp(visualizefits,'on')
        hold on
        if magictoo 
            plot(tvec{di}, X(:,ii), 'LineWidth',1,'Color',colors(1,:))
            plot(tvec{di}, X_mag(:,ii), 'LineWidth',1,'Color',colors(2,:))
            plot(t, fun(fitparam), 'LineWidth',3,'Color',colors(1,:))
            plot(t, fun(fitparam_mag), 'LineWidth',3,'Color',colors(2,:))
            legendstr = {'data','data diffused','fit','fit diffused'};
            legend(legendstr,'FontSize',20)
            titlestr = ' ';
        else
            col1faint = (colors(1,:) + 1)/(max(colors(1,:))+1);
            h3 = area(tvec{di},X(:,ii),-0.5,'FaceColor',col1faint);
            h1 = plot(tvec{di}, X(:,ii), 'LineWidth',2,'Color',colors(1,:));
            h2 = plot(t, fun(fitparam), 'LineWidth',5,'Color','k');
            set(gca,'layer','top');
            legendstr = {'data','fit','integral'};  
            legend([h1 h2 h3],legendstr,'FontSize',22)
            titlestr = ' ';%sigmoidal fit';
        end
        hold off
        xlim([min(tvec{di}),max(tvec{di})])
        ylim([-0.5,2.5])
        xlabel('time (hr)');
        yl = ylabel('SMAD4 (N:C)');
        set(yl, 'Units', 'Normalized', 'Position', [-0.18, 0.5, 0]);
        axis square;
        %title(num2str(ii))
        fp = fitparam_mag;
        I = mean(X_mag(:,ii));
        %titlestr = ['H:' num2str(fp(1),2) ', L:' num2str(fp(1)-fp(2),2) ', D:' num2str(fp(3),2) ', I:' num2str(I,2)];
        title(titlestr);
        lw = 2;
        cleanSubplot(fs,lw)
        drawnow
%         if magictoo
%             savefigure(fullfile(savedir,[dsetlabel{di}(1:6) '_fit' num2str(ii)]));
%         else
%             savefigure(fullfile(savedir,[dsetlabel{di}(1:6) '_fitRaw' num2str(ii)]));
%         end
        pause(0.7)
        
        %close
    end
end
toc
%close all

%% single-variable cell fate

chan = [1,3]; %channels to use for the log expression ratio; 1 = ISL1, 2 = SOX2, 3 = NANOG
%rthresh = [0.5, 0, 0, -0.5, -0.25]; %ratio threshold to separate ISL1/NANOG
rthresh = [0 0 0 0 0];
ratio = {};

datasetindices = 1:4;

for di = datasetindices% 4 5]

    dsetlabel{di}
    if length(chan) == 2
        %ratio = log( ( fm(chan(1),:) - min(fm(chan(1),:)) + 1) ./ ( fm(chan(2),:) - min(fm(chan(2),:)) + 1));

        % BELOW: add a number ~ 10% of minimum for regulating small values
        % makes fate distribution look much cleaner but makes no difference for
        % ordering by cell fate since I guess it is monotonic so doesn't lead
        % to reordering
        % need to exclude margin around origin
        %fm1 = fm1  - min(fm1)*0.9;
        %fm2 = fm2  - min(fm2)*0.9;
        %--------------------------------
        % removed because cannot do the same for magic
        % since that was already run on log transformed values
        %--------------------------------
        fm1 = fm{di}(chan(1),:);
        fm2 = fm{di}(chan(2),:);
        %fm1 = fm1/mean(fm1);
        %fm2 = fm2/mean(fm2);
        ratio{di} = log(fm1./fm2);

%         fm1 = fm_dmagic_fate{di}(chan(1),:);
%         fm2 = fm_dmagic_fate{di}(chan(2),:);
%         ratio{di} = fm1 - fm2;
        
        rlabel =...
            strcat("log( ",channelLabels{di}{chan(1)}," / ",channelLabels{di}{chan(2)}," )");
    else
        ratio{di} = log(fm{di}(chan,:));
        rlabel = strcat("log( ",channelLabels{di}{chan}," )");
    end

    figure,
    %figure('Position',figpos)
    [f,xi] = ksdensity(ratio{di});
    dx = xi(3)-xi(2);
    lw = 3;
    hold on
    counts = histc(ratio{di}, xi);
    bar(xi, counts/sum(counts.*dx));
    plot(xi,f,'LineWidth',4)
    hold off
    xlim([-3 3]);
    %title('Fate marker variable'); 
    xlabel(rlabel); 
    ylabel('density')
    lw = 2;
    cleanSubplot(fs,lw); axis square;
    box off
    xline(rthresh(di),'LineWidth',1.5,'LineStyle','--');
    savefigure(fullfile(savedir,[dsetlabel{di} '_fateseparation']))
    close;

%     figure('Position',figpos)
%     colorscatter(fm{di}(chan(1),:),fm{di}(chan(2),:),ratio{di},tol)
%     axlims(fm{di}(chan(1),:),fm{di}(chan(2),:),0.005);
%     cleanSubplot(fs); axis square; colormap jet
%     xlabel(channelLabels{di}{chan(1)}); ylabel(channelLabels{di}{chan(2)})
%     h = colorbar; set(get(h,'label'),'string',rlabel);
%     savefigure(fullfile(savedir,[dsetlabel{di} '_cellfatescatterplot']))
%     close;
end

%% signaling heatmap -> sorted by cell fate

idx = {};
datasetindices = 1:4;
S = {};
S_raw = {};
S_magic = {};
S_dmagic_hist = {};
S_dmagic_fate = {};

for filtertype = 1:4
    for di = datasetindices

        %------------------------------------------------------------------
        % SET IDX HERE - PAY ATTENTION
        %------------------------------------------------------------------
        idx{di} = p{di}>0;%==1; %& ~excludeidx;%p > 0; %p==4;


        if filtertype == 1
            suffix = '';
            X = hist_raw{di};
            fatem = fm{di};
            lfm = log(fatem');%./mean(fatem'));

            Xp = X(:,idx{di});
            [~, S_raw{di}, ~, ~, explained_raw{di}, ~] = pca(Xp');
            S{di} = S_raw{di};

        elseif filtertype == 2
            suffix = '_magic';
            X = hist_magic{di};
            fatem = fm{di};
            lfm = log(fatem');%./mean(fatem'));

            Xp = X(:,idx{di});
            [~, S_magic{di}, ~, ~, explained_magic{di}, ~] = pca(Xp');
            S{di} = S_magic{di};

        elseif filtertype == 3
            suffix = '_dmagic_hist';
            X = hist_dmagic_hist{di};
            lfm = fm_dmagic_hist{di}';

            Xp = X(:,idx{di});
            [~, S_dmagic_hist{di}, ~, ~, explained_dmagic_hist{di}, ~] = pca(Xp');
            S{di} = S_dmagic_hist{di};

        elseif filtertype == 4
            suffix = '_dmagic_fate';
            X = hist_dmagic_fate{di};
            lfm = fm_dmagic_fate{di}';

            Xp = X(:,idx{di});
            [~, S_dmagic_fate{di}, ~, ~, explained_dmagic_fate{di}, ~] = pca(Xp');
            S{di} = S_dmagic_fate{di};
        end

            % TICK MARKS ON AXIS ARE COMPLETELY WRONG!!
        dsetlabel{di}
        smin = 0; % 0.55
        smax = 2;

        % exclude cells near origin : has little effect
        %cutoff = 500;
        %excludeidx = (fm(chan(1),:) - min(fm(chan(1),:)) < cutoff) & (fm(chan(2),:) - min(fm(chan(2),:)) < cutoff);

        Xp = (Xp - lowSig{di})/(highSig{di}-lowSig{di});
        ratiop = lfm(idx{di},chan(1)) - lfm(idx{di},chan(2));
        [B,I] = sort(ratiop);

    %     figure('Position',figpos)
    %     [f,xi] = ksdensity(ratiop);
    %     plot(xi,f,'LineWidth',lw)
    %     xlim([-7 7]);
    %     title('Fate marker variable'); xlabel(rlabel); ylabel('density')
    %     cleanSubplot(fs); axis square
    %     xline(rthresh(di),'LineWidth',1.5,'LineStyle','--');
    %     savefigure(fullfile(savedir,[dsetlabel{di} '_fateseparation_p']))
    %     close;

        figure,%('Position',figurePosition([400 360]))
        % imagesc([min(tvec),max(tvec)],[ratio(I(1)),ratio(I(end))],X(:,I)')
        sortedX = Xp(:,I)';
        sortedX = medfilt2(sortedX,[2 2]);
        imagesc(tvec{di}, B ,sortedX)
        lw = 2; 
        fs = 26;
        colormap(flipud(brewermap([],"RdBu"))); %turbo; 
        cleanSubplot(fs,lw); 
        %axis square;
        pbaspect([0.9 1 1])
        h = colorbar('Ticks',[0 2]); set(get(h,'label'),'string','SMAD4 (N:C)');
        xlabel('time (hr)'); 
        yl = ylabel(rlabel); 
        set(yl, 'Units', 'Normalized', 'Position', [-0.17, 0.5, 0]);
        caxis([smin,smax])
        savefigure(fullfile(savedir,[dsetlabel{di} '_signalHeatmapByFate' suffix]))
        close;
    end
end
%yticks = round(linspace(1,length(ratiop),5));
%yticks = round(min(ratiop)):round(max(ratiop));

%set tick labels
%set(gca, 'YTick', yticks, 'YTickLabel', num2str(B(yticks)',2));
%set(gca, 'YTick', yticks, 'YTickLabel', string(yticks));

%scale axes so the colorbar + label is within the figure
%scaleaxes(0.975,1)

%daspect([1 0.25 1])


%% use surf to plot unequally spaced bins

figure,

% Xp = hist_raw{di}(:,idx{di});
% Xp = (Xp - lowSig)/(highSig-lowSig);
% sortedX = Xp(:,I)';
% sortedX = medfilt2(sortedX,[2 2]);
% %sortedX = imfilter(sortedX, [1 1 1]'/3);
[Xm,Ym] = meshgrid(tvec{di}, B);
surf(Xm,Ym,sortedX,'EdgeColor','none')
view([0 0 -1])
colormap(flipud(brewermap([],"RdBu"))); %turbo; 
cleanSubplot(fs,lw); 
axis square;
ylim([min(ratiop) max(ratiop)])
set(gca,'layer','top');
h = colorbar; 
set(get(h,'label'),'string','SMAD4 (N:C)');
xlabel('time (hr)'); ylabel(rlabel); 
box on
caxis([smin smax])
savefigure(fullfile(savedir,[dsetlabel{di} '_signalHeatmapByFate_uniformfateaxis']))



%% signaling heatmap -> sorted by mean signaling

s = sum(Xp,1);
[B,I] = sort(s);

figure('Position',figpos)
% imagesc([min(tvec),max(tvec)],[ratio(I(1)),ratio(I(end))],X(:,I)')
sortedX = Xp(:,I)';
%sortedX = medfilt2(sortedX,[5 3]);
imagesc([min(tvec{di}),max(tvec{di})],[min(s),max(s)],sortedX)
colormap turbo; cleanSubplot(fs); 
axis square;
h = colorbar; set(get(h,'label'),'string','SMAD4 (N:C)');
xlabel('Time (hours)'); 
ylabel('mean signaling'); caxis([smin,smax])

yticks = round(linspace(min(s),max(s),10));
%yticks = round(min(ratiop)):round(max(ratiop));
set(gca, 'YTick', yticks)

%scale axes so the colorbar + label is within the figure
%scaleaxes(0.975,1)

%daspect([1 0.25 1])

savefigure(fullfile(savedir,[dsetlabel{di} '_signalHeatmapBySignal']))

% %% signaling heatmap -> hierarchical clustering
% %%%%% TODO: double check that this is right %%%%%%%
% %methods = average, centroid, complete (farthest), median,
% %single (nearest), ward (minimum variance), weighted
% method = 'average';
% dist = 'euclidean';
% 
% Y = pdist(Xp',dist); %******X or X'?
% Z = linkage(Y,method);
% 
% figure('Position',figpos)
% % cutoff = median([Z(end-k+1,3) Z(end-k+2,3)]);
% [~, ~, outperm] = dendrogram(Z,0);%,'ColorThreshold',cutoff)
% close;
% % title('dendrogram')
% % cleanSubplot(fs); axis square
% % savefigure(fullfile(savedir,'hierarchicalDendogram'))
% 
% figure('Position',figpos)
% imagesc([min(tvec),max(tvec)],[1,length(ratiop)],Xp(:,outperm)')
% colormap turbo; cleanSubplot(fs); axis square;
% h = colorbar; set(get(h,'label'),'string','SMAD4 (N:C)');
% xlabel('Time (hours)'); caxis([smin, smax])
% 
% %scaleaxes(0.975,1)
% 
% %savefigure(fullfile(savedir,strcat([dset(1:6) '_signalHeatmapHierarchical_'],method,'_',dist)))


%% average signaling within fates

di = 3;
figure,
X = hist_raw{di}; 
X = (X - lowSig{di})/(highSig{di}-lowSig{di});
r = ratio{di};

colors = lines(2);

marg = 0.5;
lw = 3;
plot(tvec{di}, mean(X,2),'LineWidth',lw,'Color','k')
hold on
plot(tvec{di}, mean(X(:,r>0+marg),2),'LineWidth',lw,'Color',colors(1,:))
plot(tvec{di}, mean(X(:,r<0-marg),2),'LineWidth',lw,'Color',colors(2,:))
hold off
lw = 2;
fs = 26;
cleanSubplot(fs,lw)
xlim([min(tvec{di}),max(tvec{di})]); axis square
xlabel('time (hr)'); ylabel('SMAD4 (N:C)')
title(' ')
box off
legendstr = {'combined','amnion','pluripotent'};
%legendstr = {'mean','ISL1>NANOG','ISL1<NANOG'};
legend(legendstr,'Location','southwest');%'northeast') %
savefigure(fullfile(savedir,[dsetlabel{di} '_avgSignalingByFate']))

%% average signaling within fates by position subplots

ps = unique(p{di});
npos = numel(ps);

leg = {'ExE','pluripotent'}; %legend text
colors = 0.8*[1 0 0; 0 1 1]; %colors for the plots
figure('Position',figurePosition([560*npos,560]))
ax = gobjects(npos,1);

marg = 0.5;

for pidx = 1:npos
    
    legid = true(1,2);
    ax(pidx) = subplot(1,npos,pidx);
    xpos = hist_raw{di}(:,p{di}==pidx); rpos = ratio{di}(p{di}==pidx);
    disp([sum(rpos >= rthresh(di) + marg), sum(rpos < rthresh(di) - marg),...
        sum((rpos > rthresh(di) - marg) & (rpos < rthresh(di) + marg))])
    
    if sum(rpos >= rthresh(di) + marg) > 10 %don't plot average signaling in a fate with less than 10 cells
        plot(tvec{di},median(xpos(:,rpos >= rthresh(di) + marg),2),'LineWidth',lw,...
            'Color',colors(1,:))
    else
        legid(1) = false;
    end
    hold on
    if sum(rpos < rthresh(di) - marg) > 10
        plot(tvec{di},median(xpos(:,rpos < rthresh(di) - marg),2),'LineWidth',lw,...
            'Color',colors(2,:))
    else
        legid(2) = false;
    end
    
    plot(tvec{di},median(xpos(:,rpos > rthresh(di) - marg & (rpos < rthresh(di) + marg)),2),'LineWidth',lw,...
            'Color','y')
    
    hold off
%     legend(strcat(rlabel," >= ", num2str(rthresh)),...
%         strcat(rlabel," < ", num2str(rthresh)),'Location','southeast')
    legend(leg(legid),'Location','southeast')
    cleanSubplot(fs); xlim([min(tvec{di}),max(tvec{di})]); axis square
    xlabel('Time (hours)'); ylabel('SMAD4 (N:C)')
    title(conditions{di}{pidx})
end
linkaxes(ax)
ylim([0.3,1.2])

savefigure(fullfile(savedir,[dsetlabel{di} '_avgSignalingByFateByPosition']))

%% average signaling within fates by position overlayed

di = 1;
ps = unique(p{di});
npos = numel(ps);

leg = {'amnion','pluripotent'}; %legend text
%colors = 0.8*[1 0 0; 0 1 1]; %colors for the plots
colors = lines(4);
figure,%('Position',figurePosition([560*npos,560]))

marg = 0.5;

hold on
for pidx = 1:4
    
    pidx
    X = hist_dmagic_hist{di};
    lfm = fm_dmagic_hist{di}(:,p{di}==pidx)';
    %X = hist_magic{di};
    X = hist_raw{di};
	lfm = log(fm{di}(:,p{di}==pidx)');
    Xp = X(:,p{di}==pidx);
    
    lowSigp = mean(mean(Xp(1:2,:)));
    highSigp = mean(mean(Xp(6:60,:)));
    
    % normalize per condition
    Xp = (Xp - lowSigp)/(highSigp-lowSigp);
    %Xp = (Xp - lowSig{di})/(highSig{di}-lowSig{di});
    legid = true(1,2);
    
    xpos = Xp; 
    rpos = lfm(:,1)-lfm(:,3);
    
    disp([sum(rpos >= rthresh(di) + marg), sum(rpos < rthresh(di) - marg),...
         sum((rpos > rthresh(di) - marg) & (rpos < rthresh(di) + marg))])
    
    if sum(rpos >= rthresh(di) + marg) > 10 %don't plot average signaling in a fate with less than 10 cells
        plot(tvec{di},mean(xpos(:,rpos >= rthresh(di) + marg),2),'-','LineWidth',lw,...
            'Color',colors(1,:));%colors(pidx,:))
%         plot(tvec{di},cumsum(mean(xpos(:,rpos >= rthresh(di) + marg),2))/ntime,'-','LineWidth',lw,...
%             'Color',colors(pidx,:));%colors(pidx,:))
        disp(['pos mean: ' num2str(mean(mean(xpos(:,rpos >= rthresh(di) + marg),2)))])
    else
        legid(1) = false;
    end

    if sum(rpos < rthresh(di) - marg) > 10
        plot(tvec{di},mean(xpos(:,rpos < rthresh(di) - marg),2),'LineWidth',lw,...
            'Color',colors(2,:));%colors(pidx,:))
%         plot(tvec{di},cumsum(mean(xpos(:,rpos < rthresh(di) - marg),2))/ntime,'LineWidth',4,...
%             'Color',colors(pidx,:));%colors(pidx,:))
        disp(['neg mean: ' num2str(mean(mean(xpos(:,rpos < rthresh(di) - marg),2)))])
    else
        legid(2) = false;
    end
    
    %plot(tvec{di},mean(xpos,2),'-','LineWidth',1,'Color',colors(pidx,:));%colors(pidx,:))
    
%     plot(tvec{di},mean(xpos(:,rpos > rthresh(di) - marg & (rpos < rthresh(di) + marg)),2),'LineWidth',lw,...
%             'Color','y')
    
%     legend(strcat(rlabel," >= ", num2str(rthresh)),...
%         strcat(rlabel," < ", num2str(rthresh)),'Location','southeast')
    legend(leg,'Location','southwest')
    fs = 26;
    cleanSubplot(fs); xlim([min(tvec{di}),max(tvec{di})]); axis square
    xlabel('time (hr)'); ylabel('SMAD4 (N:C)')
    title(' ')
end
ylim([-0.1 , 1.8])

savefigure(fullfile(savedir,[dsetlabel{di} '_avgSignalingByFateByPosition']))


%% separate along principle components
%for previous results, either use 210801 dataset or only position/condition
%1 (30k/BMP50) for 210827 dataset

%datasetindices = 1;

S = {};
explained = {};

S_dmagic_fate = {};
explained_dmagic_fate = {};

S_dmagic_hist = {};
explained_dmagic_hist = {};

lw = 3;
fs = 26;

for di = datasetindices

    Xp = hist_raw{di}(:,idx{di});
    Xp = (Xp - lowSig{di})/(highSig{di}-lowSig{di});
    [C, S{di}, L, tsquared, explained{di}, mu] = pca(Xp');
    %X = hist{di};
    X = Xp;

    % C are the principal component vectors, column 1 is PC1, etc
    % so coefficients of PC in original basis

    % S = score are just the inner product between PC and observations?
    % so coefficients of each observation in the basis of PCs


    close all
    % figure('Position',[0 225 1500 500])
    figure('Position',figurePosition([560*3 560]))
    for ii = 1:3
        subplot(1,3,ii); hold on
        s = S{di}(:,ii);
        mu = mean(s);
        sigma = std(s);
        plot(tvec{di},mean(X(:,s < mu - sigma),2),'LineWidth',lw)
        plot(tvec{di},mean(X(:,s > mu + sigma),2),'LineWidth',lw)
        plot(tvec{di},mean(X,2),'LineWidth',lw)
        hold off
        lstr = sprintf('x_%d',ii);
        legend([lstr,'<-\sigma'],[lstr,'>\sigma'],'mean')
        title(sprintf('Principal component %d',ii))
        cleanSubplot(30); axis square %cleanSubplot(fs);
        xlim([0,42]); ylim([0 2.5])
        xlabel('time (hr)')
        if ii == 1
            ylabel('SMAD4 (N:C)')
        end
    end
    savefigure(fullfile(savedir,[dsetlabel{di} '_signalingAlongPCs']))
    close;
    
    % individual plots for paper figure
    for ii = 1:3
        figure,
        hold on
        s = S{di}(:,ii);
        mu = mean(s);
        sigma = std(s);
        plot(tvec{di},mean(X(:,s < mu - sigma),2),'LineWidth',lw)
        plot(tvec{di},mean(X(:,s > mu + sigma),2),'LineWidth',lw)
        plot(tvec{di},mean(X,2),'LineWidth',lw, 'Color','k')
        hold off
        lstr = sprintf('x_%d',ii);
        legend([lstr,'<-\sigma'],[lstr,'>\sigma'],'mean')
        %title(sprintf('PC %d (%d%%)',ii, round(explained{di}(ii))))
        %title(sprintf('principal component %d',ii))
        title(' ')
        cleanSubplot(fs); axis square %cleanSubplot(fs);
        xlim([0,42]); 
        ylim([-0.5 2.5])
        xlabel('time (hr)')
        if ii == 1
            ylabel('SMAD4 (N:C)')
        end
        savefigure(fullfile(savedir,[dsetlabel{di} '_signalingAlongPC' num2str(ii)]))
        close;
    end

    % plot histories sorted by PC

    ii = 1;
    [B,I] = sort(S{di}(:,ii),'ascend');

    figure('Position',figpos)
    % imagesc([min(tvec),max(tvec)],[ratio(I(1)),ratio(I(end))],X(:,I)')
    sortedX = Xp(:,I)';
    %sortedX = medfilt2(sortedX,[5 3]);
    imagesc([min(tvec{di}),max(tvec{di})],[min(S{di}(:,ii)),max(S{di}(:,ii))],sortedX)
    colormap turbo; cleanSubplot(fs); 
    axis square;
    h = colorbar; set(get(h,'label'),'string','SMAD4 (N:C)');
    xlabel('Time (hours)'); 
    ylabel(['PC' num2str(ii)]); caxis([smin,smax])

    %yticks = round(linspace(min(S(:,ii)),max(S(:,ii)),10));
    yticks = round(min(S{di}(:,ii))):round(max(S{di}(:,ii)));
    %yticks = round(min(ratiop)):round(max(ratiop));
    set(gca, 'YTick', yticks)
    close;
end


%% determine duration, high level, low level, integral, duration over threshold

%datasetindices = 4;%2;%[1 2 4 5];
%[1 2];
allstats = {};
ratio_raw = {};
ratio_dmagic_hist = {};
ratio_dmagic_fate = {};
ratio_dmagic_scramble = {};
confM = {};
S = {};

%%
datasetindices = 3;

for filtertype = [1 4]%[1 4]%1:4%[1 2 3 4 5]%1:3%:3

    stats = {};
    maxaccs = {}; % for classification via thresholding of signaling stats
    useValidated = false; %don't have validated histories included for now

    for di = datasetindices

        %------------------------------------------------------------------
        % SET IDX HERE - PAY ATTENTION
        %------------------------------------------------------------------
        idx{di} = p{di}>0;
        
        if filtertype == 1
            suffix = [''];
            X = hist_raw{di};
            fatem = fm{di};
            lfm = log(fatem');%./mean(fatem'));
            ratio_raw{di} = lfm(:,1) - lfm(:,3);
            ratio = ratio_raw{di};
            
            Xp = X(:,idx{di});
            [~, S_raw{di}, ~, ~, explained_raw{di}, ~] = pca(Xp');
            S{di} = S_raw{di};
            
        elseif filtertype == 2
            suffix = '_magic';
            X = hist_magic{di};
            fatem = fm{di};
            lfm = log(fatem');%./mean(fatem'));
            ratio = ratio_raw{di};
            
            Xp = X(:,idx{di});
            [~, S_magic{di}, ~, ~, explained_magic{di}, ~] = pca(Xp');
            S{di} = S_magic{di};
            
        elseif filtertype == 3
            suffix = '_dmagic_hist';
            X = hist_dmagic_hist{di};
            lfm = fm_dmagic_hist{di}';
            ratio_dmagic_hist{di} = lfm(:,1) - lfm(:,3);
            ratio = ratio_dmagic_hist{di};
            
            Xp = X(:,idx{di});
            [~, S_dmagic_hist{di}, ~, ~, explained_dmagic_hist{di}, ~] = pca(Xp');
            S{di} = S_dmagic_hist{di};
            
        elseif filtertype == 4
            suffix = ['_dmagic_fate'];
            X = hist_dmagic_fate{di};
            lfm = fm_dmagic_fate{di}';
            ratio_dmagic_fate{di} = lfm(:,1) - lfm(:,3);
            ratio = ratio_dmagic_fate{di};
            
            Xp = X(:,idx{di});
            [~, S_dmagic_fate{di}, ~, ~, explained_dmagic_fate{di}, ~] = pca(Xp');
            S{di} = S_dmagic_fate{di};
            
        elseif filtertype == 5
            suffix = '_scramble';
            X = hist_dmagic_scramble{di};
            lfm = fm_dmagic_scramble{di}';
            ratio_dmagic_scramble{di} = lfm(:,1) - lfm(:,3); 
            ratio = ratio_dmagic_scramble{di};
            
            [~, S_dmagic_scramble{di}, ~, ~, explained_dmagic_scramble{di}, ~] = pca(Xp');
            S{di} = S_dmagic_scramble{di};
        end

        % NORMALIZE
        X = (X - lowSig{di})/(highSig{di}-lowSig{di});
      
        %-----FIT------------------------------------------------------    
        
        nhists = size(X,2);
        startPoint = [0.5 1 26 1];
        Tm = zeros(nhists,1); Us = Tm; Ls = Tm; Ds = Tm; fitint = Tm;
        tic
        for ii = 1:nhists
            
            t = tvec{di}(4:end); % exclude time before BMP
            fun = @(x) x(1) - x(2)./(1 + exp(-(t-x(3))/x(4)));
            fundiff = @(x) fun(x) - X(4:end,ii);
            x0 = startPoint;
            lb = [0 0 0 0];
            ub = [2 2 42 42];
            fitparam = lsqnonlin(fundiff, x0, lb, ub);
            %[x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(___)
            Us(ii) = fitparam(1); Ls(ii) = fitparam(1)- fitparam(2); Tm(ii) = fitparam(3); 
            % if there is no difference between high/low then duration
            % is not defined 
            if Us(ii) - Ls(ii) < 0.05 
                Tm(ii) = NaN;
            end
            fitint(ii) = mean(fun(fitparam));
                
            visualizefits = 'off';
            if strcmp(visualizefits,'on')
                clf
                hold on
                plot(tvec{di}, X(:,ii), 'LineWidth',1,'Color',colors(1,:))
                %plot(tvec{di}, X_mag(:,ii), 'LineWidth',1,'Color',colors(2,:))
                plot(t, fun(fitparam), 'LineWidth',3,'Color',colors(1,:))
                %plot(t, fun(fitparam_mag), 'LineWidth',3,'Color',colors(2,:))
                hold off
                legend({'data','fit'},'FontSize',20)
                xlim([min(tvec{di}),max(tvec{di})])
                ylim([0,2.5])
                xlabel('time (hr)');
                ylabel('SMAD4 (N:C)');
                axis square;
                %title(num2str(ii))
                title(' ')
                cleanSubplot(fs)
                drawnow
                %savefigure(fullfile(savedir,[dsetlabel{di}(1:6) '_fit' num2str(ii)]));
                pause(0.7)

                %close
            end
        end
        toc

        %collect signaling features in a single matrix of size (ncells x nfeatures)
        signalthreshold = 0.75; %only consider signaling high above this level -> can make variable
        tmax = max(tvec{di});%30; % (hrs)
        rawint = mean(X(tvec{di} < tmax,:),1)'; % is nearly identical to fitint
        stats{di} = [Tm Us Ls rawint sum(X > signalthreshold,1)'*(tvec{di}(2) - tvec{di}(1))];
        stitles = {'duration (hr)',' high level','low level','SMAD4 integral','duration over threshold'};
        nstats = size(stats{di},2);
        
        %----------------------------------------------------------------
        
        allstats{filtertype, di} = stats{di};
        
        goodstatsidx = ~any(isnan(stats{di}(idx{di},:)),2);
        
        % FEATURE VS FATE
        %--------------------------------------
        
        figure('Position',[0 0 5*400 400])
        
        xlims = {[5 42],[0.7 1.2],[0.5 0.9],[0.6 1.1],[5 42]};
        for ii = 1:nstats
            subplot(1,nstats,ii); hold on
            r = ratio(idx{di});
            scatter(stats{di}(idx{di}, ii),r,100,'.')
            st = stats{di}(idx{di}, ii);
            corrcoeff = corr(st(goodstatsidx),r(goodstatsidx));
            title(corrcoeff,'FontSize',20);
            xlabel(stitles{ii})
            ylabel('log( ISL1/NANOG )')
            %xlim(xlims{ii})
            %ylim([-4 4]);
            cleanSubplot(fs); 
            axis square;
        end
        savefigure(fullfile(savedir,[dsetlabel{di} '_featureVsFate' suffix '.png']))
        close
        
        for ii = 1:nstats
            figure,
            r = ratio(idx{di});
            scatter(stats{di}(idx{di}, ii),r,100,'.')
            st = stats{di}(idx{di}, ii);
            corrcoeff = corr(st(goodstatsidx),r(goodstatsidx));
            %text(0,0.5, ['corr: ' num2str(corrcoeff,2)],'FontSize',20);
            xlabel(stitles{ii})
            ylabel('log( ISL1/NANOG )')
            %xlim(xlims{ii})
            %ylim([-4 4]);
            fs = 26;
            cleanSubplot(fs); 
            axis square;
            savefigure(fullfile(savedir,[dsetlabel{di} '_featureVsFate_' stitles{ii} '_' suffix '.png']))
            close
        end
        
        % for double magic also plot with corrected 
        if filtertype > 2
            r = lfm(:,1)-lfm(:,3);
            r = r(idx{di});
            figure('Position',[0 0 5*400 400])
            %idx = p{di}>0;
            xlims = {[5 42],[0.7 1.2],[0.5 0.9],[0.6 1.1],[5 42]};
            for ii = 1:nstats
                subplot(1,nstats,ii); hold on
                scatter(stats{di}(idx{di}, ii),r,100,'.')
                st = stats{di}(idx{di}, ii);
                corrcoeff = corr(st(goodstatsidx),r(goodstatsidx));
                title(corrcoeff,'FontSize',20);
                xlabel(stitles{ii})
                ylabel('log( ISL1/NANOG )')
                %xlim(xlims{ii})
                %ylim([-4 4]);
                cleanSubplot(fs); 
                axis square;
            end
            savefigure(fullfile(savedir,[dsetlabel{di} '_featureVsFate' suffix 'corrfate.png']))
            close
        end
        
        for ii = 4%1:4
            figure,
            hold on
            
            % also scatter raw
            lfmraw = log(fm{di}');
            rraw = lfmraw(:,1)-lfmraw(:,3);
            straw = allstats{1,di}(:,ii);
            scatter(straw, rraw,'o','MarkerEdgeColor', lines(1))
            %[corrcoeff,pval] = corr(straw(goodstatsidx),rraw(goodstatsidx));
            
            % scatter denoised
            r = lfm(:,1)-lfm(:,3);
            r = r(idx{di});
            st = stats{di}(idx{di}, ii);
            scatter(st,r,100,'.','MarkerEdgeColor', lines(1))
            [corrcoeff,pval] = corr(st(goodstatsidx),r(goodstatsidx));
            %xt = 0.95; yt = -1.35;
            xt = 0.57; yt = 1.9;
            text(xt,yt, ['corr= ' num2str(corrcoeff,2)],'FontSize',20,'HorizontalAlignment','Left');
            text(xt,yt-0.25, ['p= ' num2str(pval,2)],'FontSize',20,'HorizontalAlignment','Left');
            
            hold off
            legend({'raw','denoised'},'Location','SouthEast');%,'NorthWest')
            
            xlabel(stitles{ii})
            ylabel('log( ISL1 / NANOG )')
            %xlim(xlims{ii})
            %ylim([-1.5 1.5]);
            xlim([0.55 1.2])
            fs = 26;
            cleanSubplot(fs); 
            axis square;
            savefigure(fullfile(savedir,[dsetlabel{di} '_featureVsFate_' stitles{ii} '_' suffix 'corrfate.png']))
            %close
        end

        % PCA VS FATE
        %--------------------------------------
        
        figure('Position',[0 0 5*400 400])
        
        xlims = {[5 42],[0.7 1.2],[0.5 0.9],[0.6 1.1],[5 42]};
        for ii = 1:3 
            subplot(1,5,ii); hold on
            r = lfm(:,1)-lfm(:,3);
            r = r(idx{di});
            scatter(S{di}(idx{di}, ii),r,100,'.')
            st = S{di}(idx{di}, ii);
            corrcoeff = corr(st(goodstatsidx),r(goodstatsidx));
            title(corrcoeff,'FontSize',20);
            xlabel(['PC' num2str(ii)])
            ylabel('log( ISL1 / NANOG )')
            %xlim(xlims{ii})
            %ylim([-4 4]);
            cleanSubplot(fs); 
            axis square;
        end
        savefigure(fullfile(savedir,[dsetlabel{di} '_PCAVsFate' suffix 'corrfate.png']))
        close
    
    % classification via thresholding of signaling statistics
    %--------------------------------------
    

        %positions to use: any subset of 1:4
        poses =  unique(p{di});%[1 2 3 4];

        disp('--------------');

        %cell fate label
        label = zeros(size(ratio'));
        label(ratio > rthresh(di)) = 1; % amnion
        label(ratio <= rthresh(di)) = 2; % pluripotent

        %exclude cells with intermediate log(ISL1/NANOG) from classification:
        %highthresh = 0.5; lowthresh = 0;
        % label(ratio < 0) = 2;
        % label(ratio > 0.5) = 1;

        % mask excludes unwanted positions or cells between two thresholds 
        % (that are inbetween the fates)
        mask = ismember(p{di},poses) & ismember(label,[1,2]);
        if useValidated && true
            mask = mask & vr;
            vlabel = " validated";
        else
            vlabel = " all histories";
        end

        % close all
        figure('Position',[0 685 1500 600])
        ax = gobjects(2,nstats);
        maxMIs = NaN(nstats,2);
        maxaccs{di} = NaN(nstats,2);

        for sidx = 1:nstats % iterate over different stats
            spos = stats{di}(mask,sidx);
            lpos = label(mask)';

            npts = 100;
            vals = linspace(min(stats{di}(:,sidx)),max(stats{di}(:,sidx)),npts);
            infos = NaN(npts,2);
            accuracies = NaN(npts,2);
            classAccuracies = NaN(npts,2);
            confM{filtertype, sidx} = NaN(2,2,npts);
            for ii = 1:npts
                ltest = 2*ones(size(lpos));
                ltest(spos >= vals(ii)) = 1;
                M = confusionFromLabels(lpos,ltest); %confusion matrix from classification results
                [MI,acc] = infoFromConfusion(M); %mutual information from confusion matrix
                [MIb,accb] = infoFromConfusionBalanced(M); %MI with classes balanced
                
                confM{filtertype, sidx}(:,:,ii) = M/sum(M(:));
                infos(ii,:) = [MI, MIb];
                accuracies(ii,:) = 100*[acc, accb];
                classAccuracies(ii,1) = 100*M(1,1)/sum(M(1,:));
                classAccuracies(ii,2) = 100*M(2,2)/sum(M(2,:));
            end

            [MIbmax,I] = max(infos(:,2));
            maxMIs(sidx,:) = [MIbmax,vals(I)];
            [accmax,I1] = max(accuracies(:,1));
            [accbmax,I2] = max(accuracies(:,2));
            disp(stitles{sidx})
            disp('confusion:');
            confM{filtertype, sidx} = confM{filtertype, sidx}(:,:,I);
            confM{filtertype, sidx};
            fprintf('max accuracy = %.2g\n',accmax)
            fprintf('class balanced = %.2g\n',accbmax)
            maxaccs{di}(sidx,:) = [accbmax,vals(I2)];
        %     maxaccs(sidx,:) = [accmax,vals(I1)];

            ax(1,sidx) = subplot(2,nstats,sidx);
            plot(vals,infos(:,1),'LineWidth',3); hold on
            plot(vals,infos(:,2),'LineWidth',3)
            xline(vals(I),'--','LineWidth',2); hold off
            ylabel('MI (bits)'); xlabel(stitles{sidx})
            xlim(vals([1,end]))
            cleanSubplot(14); axis square
            if sidx == nstats
                legend('raw','class balanced','Location','southeast')
            end

            ax(2,sidx) = subplot(2,nstats,sidx + nstats);
            plot(vals,classAccuracies(:,1),'LineWidth',3); hold on
            plot(vals,classAccuracies(:,2),'LineWidth',3)
            xline(vals(I),'--','LineWidth',2); hold off
            ylabel('accuracy (%)'); xlabel(stitles{sidx})
            cleanSubplot(14); axis square
            ylim([0,100])
            xlim(vals([1,end]))
            if sidx == nstats
                legend('ExE','pluripotent','Location','southwest')
            end

        %     ax(2,sidx) = subplot(2,nstats,sidx + nstats);
        %     plot(vals,accuracies(:,1),'LineWidth',3); hold on
        %     plot(vals,accuracies(:,2),'LineWidth',3)
        %     xline(vals(I),'--','LineWidth',2); hold off
        %     ylabel('accuracy (%)'); xlabel(stitles{sidx})
        %     cleanSubplot; axis square
        %     ylim([0,100])

        end
        linkprop(ax(1,:),'Ylim');
        condstring = cell2mat(strcat(conditions{di}(poses),"; "));
        sgtitle(strcat(condstring, vlabel),...
            'FontSize',fs,'FontWeight','bold')
        disp(maxaccs{di})

        savefigure(fullfile(savedir,[dsetlabel{di} '_classificationBySignalFeatures' suffix]))
        close;
    
    % DISTRIBUTIONS
    %--------------------------------------
        
        ms = 40; % marker size for data scater
        
        for conditional = [false, true]
            
            % make and display joint probability plots
            nstats = 4;%size(stats{di},2);
            % figure('Position',[0 250 1500 375])
            separatedists = true;
            if ~separatedists 
                figure('Position',figurePosition([(nstats)*560,560]))
            end
            for sidx = 1:nstats
                if separatedists
                    figure,
                    %     figure('Position',figpos)
                else
                    %subplots of one figure
                    subplot(1,nstats,sidx)
                end

                np = 100;
    %             ydata = ratio{di}(:); 
    %             xdata = stats{di}(:,sidx);

                ydata = lfm(:,chan(1))-lfm(:,chan(2)); 
                xdata = stats{di}(:,sidx);

                [xl,yl] = axlims(xdata,ydata,0.02,false); % 0.0
                [XX,YY] = meshgrid(linspace(xl(1),xl(2),np),linspace(yl(1),yl(2),np));
                [f,xi] = ksdensity([xdata,ydata],[XX(:),YY(:)]);

                F = reshape(f,[np,np]);
                xx = reshape(xi(:,1),[np,np]);
                yy = reshape(xi(:,2),[np,np]);
                F = F./sum(F(:));
                F2 = F./sum(F, 1);
                F2 = F2./sum(F2(:));

                t = maxaccs{di}(sidx,2); % feature threshold
                ft = 0;% fate threshold

                truepos = sum(sum(F.*(xx>maxaccs{di}(sidx,2)).*(yy>ft)));
                falsepos = sum(sum(F.*(xx>maxaccs{di}(sidx,2)).*(yy<ft)));
                trueneg = sum(sum(F.*(xx<maxaccs{di}(sidx,2)).*(yy<ft)));
                falseneg = sum(sum(F.*(xx<maxaccs{di}(sidx,2)).*(yy>ft)));
                M = [falseneg, truepos; trueneg, falsepos];

                disp('----------')
                truepos = sum(sum(F2.*(xx>maxaccs{di}(sidx,2)).*(yy>ft)));
                falsepos = sum(sum(F2.*(xx>maxaccs{di}(sidx,2)).*(yy<ft)));
                trueneg = sum(sum(F2.*(xx<maxaccs{di}(sidx,2)).*(yy<ft)));
                falseneg = sum(sum(F2.*(xx<maxaccs{di}(sidx,2)).*(yy>ft)));
                M2 = [falseneg, truepos; trueneg, falsepos];

                if conditional
                    disttype = 'conditional';
                    F = F2;
                    M = M2;
                else
                    disttype = 'joint';
                end
            
                % override with raw data confusion matrix
                M = confM{filtertype, sidx};
                M = M(:,[2 1]);
                
                surf(xx,yy,F,'LineStyle','none')
                hold on 
                scatter3(xdata, ydata, xdata*0+ max(F(:)),60,'.','MarkerEdgeColor','#777777');
                
                % plot raw on top
                ydata = ratio_raw{di};
                xdata = allstats{1, di}(:,sidx);
                scatter3(xdata, ydata, xdata*0+ max(F(:)),30,'o','MarkerEdgeColor','#777777');

                fs = 26;
                view(2); 
                cleanSubplot(fs); 
                axis square; 
                %colormap jet
                colormap(flipud(brewermap([],"RdBu")))
                %xlims = {[12 26],[0.85 1.05],[0.65 0.75],[0.7 0.9],[15 35]};
                %xlim(xlims{sidx}); ylim([-3 3])
                xlim(xl); ylim(yl)
                if (sidx == 1 && ~separatedists) || (sidx == 3 && separatedists)
                    %ylabel(strcat("fate = ",rlabel))
                    ylabel(rlabel)
                end
                set(gca, 'YTick', [-1 0 1]);
                xlabel(stitles{sidx})

                lw = 3;
                xline(t,'--','LineWidth',lw,'Color','k');
                yline(ft,'--','LineWidth',lw);

                title(' ');
                set(gca,'layer','top');
                grid off
                box on

                if separatedists
                    savefigure(fullfile(savedir,[dsetlabel{di} '_signalFeaturesDistributions_' disttype '_' stitles{sidx} suffix]))
                end

                zc = max(F(:))+0.01;
                zlim([0 zc]);

                fs = 22;
                cb ='k';
                 text(xl(1) + (t-xl(1))/2, yl(2)/2, zc, [num2str(round(M(1,1)*100)) '%'], 'Color',cb,'FontSize',fs,'HorizontalAlignment','center','FontWeight','bold');%,'BackgroundColor','white')
                text(t + (xl(2)-t)/2, yl(2)/2, zc, [num2str(round(M(1,2)*100)) '%'], 'Color',cb,'FontSize',fs,'HorizontalAlignment','center','FontWeight','bold')
                text(xl(1) + (t-xl(1))/2, yl(1)/2, zc, [num2str(round(M(2,1)*100)) '%'], 'Color',cb,'FontSize',fs,'HorizontalAlignment','center','FontWeight','bold')
                text(t + (xl(2)-t)/2, yl(1)/2, zc, [num2str(round(M(2,2)*100)) '%'], 'Color',cb,'FontSize',fs,'HorizontalAlignment','center','FontWeight','bold')

                fs = 20;
                cnb = 'w';
                text(xl(1) + (t-xl(1))/2, yl(2)/2, zc, [num2str(round(M(1,1)*100)) '%'], 'Color',cnb,'FontSize',fs,'HorizontalAlignment','center');
                text(t + (xl(2)-t)/2, yl(2)/2, zc, [num2str(round(M(1,2)*100)) '%'], 'Color',cnb,'FontSize',fs,'HorizontalAlignment','center');
                text(xl(1) + (t-xl(1))/2, yl(1)/2, zc, [num2str(round(M(2,1)*100)) '%'], 'Color',cnb,'FontSize',fs,'HorizontalAlignment','center');
                text(t + (xl(2)-t)/2, yl(1)/2, zc, [num2str(round(M(2,2)*100)) '%'], 'Color',cnb,'FontSize',fs,'HorizontalAlignment','center');

                title(sprintf('accuracy = %.2g%%; MI = %.2g',maxaccs{di}(sidx,1),maxMIs(sidx,1)))
                hold off

                %fs = 26;
                %cleanSubplot(fs,lw)

                if separatedists
                    savefigure(fullfile(savedir,[dsetlabel{di} '_signalFeaturesDistributions_' disttype '_annot_' stitles{sidx} suffix]))
                    close;
                end
            end
            if ~separatedists
                savefigure(fullfile(savedir,[dsetlabel{di} '_signalFeaturesDistributions' disttype suffix]))
            end
            close;
        end
        
        % SCATTER OF PCs Vs FEATURES
        %--------------------------------------

        % % correlation between PCs is impossible by definition
        % scatter(S(:,1), S(:,2))
        % axis square;

        % but do our fitted coefficients correlate with PCs?
        figure('Position',figurePosition([400*3 400*4]))
        clf
        for ii = 1:4
            for jj = 1:3 % PCs
                subplot(4,3, 3*(ii-1) + jj); 
                st = stats{di}(idx{di},ii);
                scatter(S{di}(:,jj),st,'.')
                ylabel(stitles{ii})
                xlabel(['PC' num2str(jj)])
                corrcoeff = corr(st(goodstatsidx), S{di}(goodstatsidx,jj));
                title(corrcoeff,'FontSize',20)
                axis square;
            end
        end
        savefigure(fullfile(savedir,[dsetlabel{di} '_featuresVsPCA' suffix]))
        close;

        %--------------------------------------
        % HEATMAP stats vs PC
        %--------------------------------------
        figpos(3:4) = [640 560];
        figure('Position',figpos)
        shortertitles = {'duration','high','low','integral'};
        if filtertype == 1
            varidx = 1:3;
        else
            varidx = 1:4;
        end
        st = stats{di}(idx{di},varidx);
        C = corr(st(goodstatsidx,:), S{di}(goodstatsidx,1:3));
        C = round(C,2);
        writematrix(C,fullfile(savedir,['featuresVsPCA_' suffix '.csv']))

        imagesc(C,[-1 1]);
        fs = 32;
        cleanSubplot(fs); 
        colormap(flipud(brewermap([],"RdBu")))
        set(gca, 'YTick', 1:size(C,1), 'YTickLabel', shortertitles);
        set(gca, 'XTick', 1:size(C,2), 'XTickLabel', {'PC1','PC2','PC3'});
        %xtickangle(gca,40)
        colorbar('Ticks',[-1 0 1]);
        title(' ')
        CasStr = arrayfun(@num2str, C, 'UniformOutput', false);
        [X,Y] = meshgrid(1:size(C,2),1:size(C,1));
        text(X(:),Y(:),CasStr(:),'FontSize',30,'HorizontalAlignment','center')

        savefigure(fullfile(savedir,[dsetlabel{di} '_featuresVsPCAHeat' suffix]))
        close;

        %  do fit features correlate with each other?
        % mild negative correlation between high and duration, low and duration
        % higher correlation between integral and high, low, than integral and
        % duration??!

        % after magic everything correlates strongly with integral 
        % weak positive correlation between high/low and duration

        % features vs features
        figure('Position',figurePosition([400*4 400*4])),
        clf
        for ii = 1:4
            for jj = ii:4 % PCs
                subplot(4,4, 4*(ii-1) + jj); 
                scatter(stats{di}(:,jj),stats{di}(:,ii),'.')
                ylabel(stitles{ii})
                xlabel(stitles{jj})
                st = stats{di}(:,ii);
                corrcoeff = corr(st(goodstatsidx), stats{di}(goodstatsidx,jj));
                title(corrcoeff,'FontSize',20)
                axis square;
            end
        end
        savefigure(fullfile(savedir,[dsetlabel{di} '_featuresVsFeatures' suffix]))
        close;

        % heatmap stats vs stats
        figpos(3:4) = [600 560];
        figure('Position',figpos)
        shortertitles = {'duration','high','low','integral'};
        varidx = 1:4;
        st = stats{di}(:,varidx);
        C = corr(st(goodstatsidx,:), st(goodstatsidx,:));
        C = round(C,2);
        imagesc(C,[-1 1]);
        cleanSubplot(fs); colormap(flipud(brewermap([],"RdBu")))
        set(gca, 'YTick', 1:size(C,1), 'YTickLabel', shortertitles);
        set(gca, 'XTick', 1:size(C,2), 'XTickLabel', shortertitles);
        xtickangle(gca,45)
        colorbar('Ticks',[-1 0 1]);
        title(' ')
        CasStr = arrayfun(@num2str, C, 'UniformOutput', false);
        [X,Y] = meshgrid(1:4,1:4);
        text(X(:),Y(:),CasStr(:),'FontSize',30,'HorizontalAlignment','center')
        savefigure(fullfile(savedir,[dsetlabel{di} '_featuresVsFeaturesHeat' suffix]))
        close;

        % heatmap stats vs fate
        figpos(3:4) = [620 560];
        figure('Position',figpos)
        varidx = 1:4;
        st = stats{di}(:,varidx);
        C = corr(st(goodstatsidx,:), lfm(goodstatsidx,:));
        C = round(C,2);

        if di == 3 && filtertype > 2
            clabels = channelLabels{di}([1:3 5 7 12 13]);
        else
            clabels = channelLabels{di};
        end
        heatmap(clabels, shortertitles,C);%,'CellLabelColor', 'None')
        colorbar off;
        set(gca,'FontSize', fs)
        savefigure(fullfile(savedir,[dsetlabel{di} '_featuresVsfateHeat' suffix]))
        close;
    end
end

%%

%--------------------------------------
% classification via thresholding of signaling statistics
%--------------------------------------

maxaccs = {};
filtertypes = {'raw','magic','dmagic_hist','dmagic_fate'};

di = 1;

for filtertype = 1:4
    
    if filtertype == 1
        ratio = ratio_raw{di};

    elseif filtertype == 2
        ratio = ratio_raw{di};

    elseif filtertype == 3
        ratio = ratio_dmagic_hist{di};

    elseif filtertype == 4
        ratio = ratio_dmagic_fate{di};
    end

    %positions to use: any subset of 1:4
    poses =  unique(p{di});%[1 2 3 4];

    disp('--------------');
    disp(filtertypes{filtertype})
    disp('--------------');

    %cell fate label
    label = zeros(size(ratio'));
    marg = 0;
    label(ratio > rthresh(di) + marg) = 1; % amnion
    label(ratio < rthresh(di) - marg) = 2; % pluripotent

    % mask excludes unwanted positions or cells between two thresholds 
    % (that are inbetween the fates)
    mask = ismember(p{di},poses) & ismember(label,[1,2]);
    if useValidated && true
        mask = mask & vr;
        vlabel = " validated";
    else
        vlabel = " all histories";
    end

    maxMIs = NaN(nstats,2);
    maxaccs{di} = NaN(nstats,2);

    for sidx = 1:nstats % iterate over different stats
        spos = allstats{filtertype, di}(mask,sidx);
        lpos = label(mask)';

        npts = 1000;
        vals = linspace(min(allstats{filtertype, di}(:,sidx)),max(allstats{filtertype, di}(:,sidx)),npts);
        infos = NaN(npts,2);
        accuracies = NaN(npts,2);
        classAccuracies = NaN(npts,2);
        for ii = 1:npts
            ltest = 2*ones(size(lpos));
            ltest(spos >= vals(ii)) = 1;
            M = confusionFromLabels(lpos,ltest); %confusion matrix from classification results
            [MI,acc] = infoFromConfusion(M); %mutual information from confusion matrix
            [MIb,accb] = infoFromConfusionBalanced(M); %MI with classes balanced

            infos(ii,:) = [MI, MIb];
            accuracies(ii,:) = 100*[acc, accb];
            classAccuracies(ii,1) = 100*M(1,1)/sum(M(1,:));
            classAccuracies(ii,2) = 100*M(2,2)/sum(M(2,:));
        end

        [MIbmax,I] = max(infos(:,2));
        maxMIs(sidx,:) = [MIbmax,vals(I)];
        [accmax,I1] = max(accuracies(:,1));
        [accbmax,I2] = max(accuracies(:,2));
        disp(stitles{sidx})
        fprintf('max accuracy = %.2g\n',accmax)
        fprintf('class balanced = %.2g\n',accbmax)
        maxaccs{di}(sidx,:) = [accbmax,vals(I2)];
    end
    disp(maxaccs{di})
end


%% scatter PC vs fitparam w/ and w/o MAGIC (fig for 220928)

figure, 
for ii = [1 4] % which fit param
    ylims = {[7 42],[0 2],[0 2],[0.4 1.75]};
    scsym = {'.','.'};
    pcai = 1; % which PC
    corrcoeff = [0 0];
    clf
    hold on
    for filtertype = [1 2]
        scatter(S_raw{di}(:,pcai),allstats{filtertype,di}(:,ii), 500,scsym{filtertype})
        x = S_raw{di}(:,pcai);
        y = allstats{filtertype,di}(:,ii);
        good = ~isnan(x) & ~isnan(y);
        corrcoeff(filtertype) = corr(x(good),y(good));
    end
    hold off
    yl = ylabel(stitles{ii});
    set(yl, 'Units', 'Normalized', 'Position', [-0.17, 0.5, 0]);
    xlabel(['principal component ' num2str(pcai)])
    ylim(ylims{ii})
    axis square;
    lw = 2;
    cleanSubplot(fs,lw)
    %legend({['raw corr: ' num2str(corrcoeff(1),2)],['denoised: ' num2str(corrcoeff(2),2)]},'Location','SouthEast');
    legend({'raw','denoised'},'Location','SouthEast');
    savefigure(fullfile(savedir,[dsetlabel{di} '_' stitles{ii} '_PC' num2str(pcai) '_rawVmagic']))
end


%% fitparam vs fate (SI fig for 220928)

ncols = 4+1;
nrows = 4;
L = 180;
L2 = L+50;
fs = 14;
figure('Position',figurePosition([L*ncols L2*nrows]))
xlims = {[10 35],[0.5 1.7],[0.2 1],[0.4 1.25]};

    
for ii = 1:4
        
    for param = 1:4

        filtertypes = [1 3 4 5];
        filtertype = filtertypes(ii);
        
        if filtertype <= 2
            r = ratio_raw{di}';
        elseif filtertype == 3
            r = ratio_dmagic_hist{di};
        elseif filtertype == 4
            r = ratio_dmagic_fate{di};
        end

        pli = (ii-1)*ncols + param;
        subplot(nrows,ncols,pli)
        
        st = allstats{filtertype,di}(:,param);
        good = ~isnan(st(:)) & ~isnan(r(:));
        [corrcoeff, pval] = corr(st(good),r(good));
        scatter(st, r, '.')
        title(['corr=' num2str(corrcoeff,2) '\newlinep=' num2str(pval,2)]);
        
        if ii < nrows
            set(gca,'XTickLabel','')
        else
            xlabel(stitles{param})
        end
        if param > 1
            set(gca,'YTickLabel','')
        else
            ylabel('log( ISL1 / NANOG )');
        end
        
        
        ylim([-2 2])
        xlim(xlims{param});
        cleanSubplot(fs,lw)
        axis square
    end
    
    pli = (ii-1)*ncols + 5;
    subplot(nrows,ncols,pli)
    hist(r,linspace(-2,2,50))
    h = findobj(gca,'Type','patch');
    h.FaceColor = lines(1);
    h.EdgeColor = 'w';
    xlabel('log( ISL1 / NANOG )');
    ylabel('cell number');
    cleanSubplot(fs,lw)
    axis square
end
savefigure(fullfile(savedir,[dsetlabel{di} '_featuresVsFateAllScatters']))


%% fitparam vs other genes

ncols = 4+1;
nrows = 4;
L = 180;
L2 = L+50;
fs = 14;
figure('Position',figurePosition([L*ncols L2*nrows]))
xlims = {[10 35],[0.5 1.7],[0.2 1],[0.4 1.25]};

    
for ii = 1:nrows
        
    for param = 1:4

        filtertypes = [1 2 3 4];
        filtertype = filtertypes(ii);
        
        if filtertype <= 2
            fatem = fm{di};
            lfm = log(fatem');
            ci = 7;
            r = lfm(:,ci);
        elseif filtertype == 3
            lfm = fm_dmagic_hist{di}';
            r = lfm(:,5);
        elseif filtertype == 4
            lfm = fm_dmagic_fate{di}';
            r = lfm(:,5);
        end

        pli = (ii-1)*ncols + param;
        subplot(nrows,ncols,pli)
        
        st = allstats{filtertype,di}(:,param);
        good = ~isnan(st(:)) & ~isnan(r(:));
        [corrcoeff, pval] = corr(st(good),r(good));
        scatter(st, r, '.')
        title(['corr=' num2str(corrcoeff,2) '\newlinep=' num2str(pval,2)]);
        
        if ii < nrows
            set(gca,'XTickLabel','')
        else
            xlabel(stitles{param})
        end
        if param > 1
            set(gca,'YTickLabel','')
        else
            ylabel(['log( ' channelLabels{di}{ci} ' )']);
        end
        
        
        %ylim([-2 2])
        %xlim(xlims{param});
        cleanSubplot(fs,lw)
        axis square
    end
    
    pli = (ii-1)*ncols + 5;
    subplot(nrows,ncols,pli)
    hist(r,50)
    h = findobj(gca,'Type','patch');
    h.FaceColor = lines(1);
    h.EdgeColor = 'w';
    xlabel(['log( ' channelLabels{di}{ci} ' )']);
    ylabel('cell number');
    cleanSubplot(fs,lw)
    axis square
end

savefigure(fullfile(savedir,[dsetlabel{di} '_featuresVs' channelLabels{di}{ci} 'AllScatters']))

%% integral per condition violin plot (fig for 210827)

di = 1;
marg = 0.5;

for filtertype = 1%:4
    
meanstat = zeros([2 numel(unique(p{di}))]);
for ii = 4

    figure, 
    if filtertype == 1
        suffix = '';
        X = hist_raw{di};
        fatem = fm{di};
        lfm = log(fatem');%./mean(fatem')); 
    elseif filtertype == 2
        suffix = '_magic';
        X = hist_magic{di};
        fatem = fm{di};
        lfm = log(fatem');%./mean(fatem'));          
    elseif filtertype == 3
        suffix = '_dmagic_hist';
        X = hist_dmagic_hist{di};
        lfm = fm_dmagic_hist{di}';     
    elseif filtertype == 4
        suffix = '_dmagic_fate';
        X = hist_dmagic_fate{di};
        lfm = fm_dmagic_fate{di}';
    end
    fate = lfm(:,1) - lfm(:,3);
    
    % normalize per condition
	for pidx = unique(p{di})
        Xp = X(:,p{di}==pidx);
        lowSigp = mean(mean(Xp(1:2,:)));
        highSigp = mean(mean(Xp(6:60,:)));
        X(:,p{di}==pidx) = (X(:,p{di}==pidx) - lowSigp)/(highSigp-lowSigp);
    end
    
    if ii == 4 % integral
        st = mean(X,1);
    else
        st = allstats{filtertype,di}(:,ii);
    end
    
    pidx = arrayfun(@num2str, p{di} + 4*(fate'>0+marg), 'UniformOutput', false);
    excludeidx = fate'<0+marg & fate'>0-marg;
    st = st(~excludeidx);
    fate = fate(~excludeidx);
    
    colors = lines(4);
    colors = cat(1, colors([2 1 3],:), colors([2 1 3 4],:));
    violinplot(st,pidx(~excludeidx),'GroupOrder',{'2','1','3','4','6','5','7','8'},'ViolinColor',colors);
%     colors = colors([2 2 1 1 3 3 4 4],:);
%     violinplot(st,pidx(~excludeidx),'GroupOrder',{'2','6','1','5','3','7','8'},'ViolinColor',colors);
    
    fs = 26;
    cleanSubplot(fs,lw)
    title(suffix)
    axis square;
    box off;
    ylim([0.25 1.6])
    yline(0.75)
    ylabel('SMAD4 integral')
   
%     for ci = unique(p{di})
%         st = allstats{filtertype,di}(:,ii);
%         meanstat(1, ci) = mean(st(p{di} == ci & fate' > 0));
%         meanstat(2, ci) = mean(st(p{di} == ci & fate' < 0));
%     end
    savefigure(fullfile(savedir,[dsetlabel{di} '_' stitles{ii} '_violin']))
end
end


%% scatter feature vs fate colored for conditions

colorbycondition = false;
di = 3;
figure, 
for ii = 1:4%[1 4] % which fit param
    ylims = {[10 35],[0.6 1.6],[0.2 1.1],[0.3 1.35]};
    scsym = {'o','.','.','.'};
    scsiz = {50,100,100,100};
    pcai = 1; % which PC
    corrcoeff = [0 0];
    clf
    hold on
    for filtertype = [1 4]%[1 3]%[1 4]
        
        if filtertype == 1
            suffix = '';
            X = hist_raw{di};
            fatem = fm{di};
            lfm = log(fatem');%./mean(fatem'));          
        elseif filtertype == 2
            suffix = '_magic';
            X = hist_magic{di};
            fatem = fm{di};
            lfm = log(fatem');%./mean(fatem'));          
        elseif filtertype == 3
            suffix = '_dmagic_hist';
            X = hist_dmagic_hist{di};
            lfm = fm_dmagic_hist{di}';     
        elseif filtertype == 4
            suffix = '_dmagic_fate';
            X = hist_dmagic_fate{di};
            lfm = fm_dmagic_fate{di}';
        end
        
        y = lfm(:,1) - lfm(:,3);%S_raw{di}(:,pcai);
        x = allstats{filtertype,di}(:,ii);
        good = ~isnan(x) & ~isnan(y);
        colormap lines
        if colorbycondition
            scatter(x,y, scsiz{filtertype}, p{di}, scsym{filtertype})
            suffix = ['condition_' suffix];
        else
            scatter(x,y, scsiz{filtertype}, lines(1), scsym{filtertype})
        end
        [corrcoeff(filtertype),pval] = corr(x(good),y(good));
        if filtertype>1
            text(ylims{ii}(2), -1.5, ['corr=' num2str(corrcoeff(filtertype),2)],'FontSize',20,'HorizontalAlignment','right');
            text(ylims{ii}(2), -1.75, ['p=' num2str(pval,2)],'FontSize',20,'HorizontalAlignment','right');
        end
    end
    hold off
    yl = ylabel('log( ISL1 / NANOG )');
    set(yl, 'Units', 'Normalized', 'Position', [-0.15, 0.5, 0]);
    xlabel(stitles{ii})
    xlim(ylims{ii})
    ylim([-2 2])
    axis square;
    lw = 2;
    fs =26;
    cleanSubplot(fs,lw)
    %legend({['raw corr: ' num2str(corrcoeff(1),2)],['denoised: ' num2str(corrcoeff(2),2)]},'Location','SouthEast');
    legend({'raw','denoised'},'Location','NorthWest');%,'SouthEast');
    savefigure(fullfile(savedir,[dsetlabel{di} '_' stitles{ii} '_fate_' suffix]))
end


%% variance explained with and without magic

colors = lines(2);
nPC = 10;
figure, 
hold on
plot(1:nPC, explained_raw{di}(1:nPC), '-','Color',colors(1,:),'LineWidth',3)
plot(1:nPC, explained_magic{di}(1:nPC), '-','Color',colors(2,:),'LineWidth',3)
scatter(1:nPC, explained_raw{di}(1:nPC), 1000, '.','MarkerEdgeColor',colors(1,:))
scatter(1:nPC, explained_magic{di}(1:nPC), 1000, '.','MarkerEdgeColor',colors(2,:))
legend({'raw','denoised'})
%scatter(1:nPC, explained_dmagic_fate{di}(1:nPC), 1000, '.')
xlabel('principal component')
yl = ylabel('var explained (%)');
title(' ')
axis square
hold off
lw = 2;
cleanSubplot(fs,lw)
set(yl, 'Units', 'Normalized', 'Position', [-0.17, 0.5, 0]);
%savefigure(fullfile(savedir,[dsetlabel{di} '_varianceExplained']))








%% ----------------------------------------------------------------------
% OLD 
% ----------------------------------------------------------------------

%% scatter PCA colored for fate discrete
figure,
di = 5;

cidx = (ratio{di}(idx{di})>rthresh(di)) + 1;
colors = [0 1 1; 1 0 0];
c = colors(cidx,:);

S = S_raw;
scatter(S{di}(:,1), S{di}(:,2),100,c,'.')
%scatter(S_dmagic_fate{di}(:,1), S_dmagic_fate{di}(:,2),70,mat2gray(ratio{di}(idx{di})>rthresh(di)),'.')
%scatter(S_dmagic_hist{di}(:,1), S_dmagic_hist{di}(:,2),70,mat2gray(ratio{di}(idx{di})>rthresh(di)),'.')
xlabel('PC1')
ylabel('PC2')
%colormap prism
cleanSubplot(30); 
xlim([-4 6])
axis square
%savefigure(fullfile(savedir,[dsetlabel{di} '_PC1_PC2_fate' suffix]))

%% scatter PCA vs fate 

figure,
di = 5;

cidx = (ratio{di}(idx{di})>rthresh(di)) + 1;
colors = [0 1 1; 1 0 0];
c = colors(cidx,:);

scatter(S_dmagic_hist{di}(:,1), ratio{di}(idx{di}),100,'.')
corr(S_dmagic_hist{di}(:,1), ratio{di}(idx{di})')
%scatter(S_dmagic_fate{di}(:,1), S_dmagic_fate{di}(:,2),70,mat2gray(ratio{di}(idx{di})>rthresh(di)),'.')
%scatter(S_dmagic_hist{di}(:,1), S_dmagic_hist{di}(:,2),70,mat2gray(ratio{di}(idx{di})>rthresh(di)),'.')
xlabel('PC1')
ylabel('log( ISL1 / NANOG )')
%colormap prism
cleanSubplot(30); 
%xlim([-4 6])
axis square

%% scatter PCA colored for fate continuous

figure,
di = 5;

cidx = (ratio{di}(idx{di})>rthresh(di)) + 1;
colors = [0 1 1; 1 0 0];
c = colors(cidx,:);

scatter(S_dmagic_hist{di}(:,1), S_dmagic_hist{di}(:,2),100,ratio{di}(idx{di}),'.')
%scatter(S_dmagic_fate{di}(:,1), S_dmagic_fate{di}(:,2),70,mat2gray(ratio{di}(idx{di})>rthresh(di)),'.')
%scatter(S_dmagic_hist{di}(:,1), S_dmagic_hist{di}(:,2),70,mat2gray(ratio{di}(idx{di})>rthresh(di)),'.')
xlabel('PC1')
ylabel('PC2')
c = colorbar;   
c.Label.String = 'fate = log( ISL1 / NANOG )';
%colormap prism
cleanSubplot(30); 
xlim([-4 6])
axis square
%savefigure(fullfile(savedir,[dsetlabel{di} '_PC1_PC2_fate' suffix]))






%% test if correlation plots match python
% THERE IS A MISTAKE, CORRELATIONS ABOVE ARE NOT THE SAME AS IN PYTHON
% SCATTER PLOTS BELOW ARE GOOD THOUGH

di = 1;
fate = fm_dmagic_fate{di}(1,:) - fm_dmagic_fate{di}(3,:); % log isl1/nanog
sig = mean(hist_dmagic_fate{di},1);

figure, scatter(sig, fate,'.') % this is showing the original fate markers, not corrected
%%
fate = fm_dmagic_hist{di}(1,:) - fm_dmagic_hist{di}(3,:); % log isl1/nanog
sig = mean(hist_dmagic_hist{di},1);

figure, scatter(sig, fate,'.') % this is showing the original fate markers, not corrected

%%

figure,
r = lfm(:,1)-lfm(:,2);
scatter(stats{di}(idx, 4),r,'.')

%lfm = fm_dmagic_fate{di}';

% it is the cell fate variable that is different
%scatter(stats{di}(idx, 4),mean(hist_dmagic_fate{di},1),'.')

            
%%
[A,B,r,U,V] = canoncorr(stats{di}(:,varidx), lfm);

t = tiledlayout(2,2);
title(t,'Canonical Scores of X vs Canonical Scores of Y')
xlabel(t,'Canonical Variables of X')
ylabel(t,'Canonical Variables of Y')
t.TileSpacing = 'compact';

nexttile
plot(U(:,1),V(:,1),'.')
xlabel('u1')
ylabel('v1')

nexttile
plot(U(:,2),V(:,1),'.')
xlabel('u2')
ylabel('v1')

nexttile
plot(U(:,1),V(:,2),'.')
xlabel('u1')
ylabel('v2')

nexttile
plot(U(:,2),V(:,2),'.')
xlabel('u2')
ylabel('v2')

%% sort tracks by fit parameters

ii = 4;
[B,I] = sort(stats{di}(:,ii),'ascend');

figure('Position',figpos)
% imagesc([min(tvec),max(tvec)],[ratio(I(1)),ratio(I(end))],X(:,I)')
sortedX = Xp(:,I)';
%sortedX = medfilt2(sortedX,[5 3]);
imagesc([min(tvec{di}),max(tvec{di})],[min(stats{di}(:,ii)),max(stats{di}(:,ii))],sortedX)
colormap turbo; cleanSubplot(fs); 
axis square;
h = colorbar; set(get(h,'label'),'string','SMAD4 (N:C)');
xlabel('Time (hours)'); 
ylabel(stitles{ii}); caxis([smin,smax])

ytickss = linspace(min(stats{di}(:,ii)),max(stats{di}(:,ii)),10);
%yticks = round(min(ratiop)):round(max(ratiop));
set(gca, 'YTick', ytickss)


%% local functions

function scaleaxes(xscale,yscale)

pos = get(gca,'Position');
pos(2) = pos(2) + (1 - yscale)*pos(4);
pos(4) = yscale*pos(4);
pos(1) = pos(1) - (1 - xscale)*pos(3);
pos(3) = xscale*pos(3);
set(gca,'Position',pos)

end

