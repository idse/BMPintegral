clear; close all; clc

%% setup
%graphical options
fs = 26; %font size
wh = [560 560]; figpos = figurePosition(wh); %figure width and height -> figure position
lw = 3; %line width

scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
dataDir = fullfile(scriptPath,'data');
savedir = fullfile(scriptPath,'figures');

%dsetlabel = {'210801_histories','210827_histories','210827_corrected', '220928_histories'};%,'210827_corrected_only_updated','211115_histories','220111_histories'};
dsetlabel = {'210827_histories_nomedfilt','220928_histories_nomedfilt'};
ndsets = length(dsetlabel);

hist_raw = cell(1,ndsets);
hist_magic = cell(1,ndsets);
% hist_dmagic = {};
hist_dmagic_hist = cell(1,ndsets);
hist_dmagic_fate = cell(1,ndsets);
hist_dmagic_scramble = cell(1,ndsets);
tvec = cell(1,ndsets);
fm = cell(1,ndsets);
% fm_dmagic = {};
fm_dmagic_scramble = cell(1,ndsets);
fm_dmagic_hist = cell(1,ndsets);
fm_dmagic_fate = cell(1,ndsets);
p = cell(1,ndsets);
conditions = cell(1,ndsets);
channelLabels = cell(1,ndsets);

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
    
    %get a vector of imaging time points for each dataset
    ntime = s.liveMeta.nTime;
    timeres = strsplit(s.liveMeta.timeInterval,' ');
    timeres = str2double(timeres{1});
    tscale = timeres/60;
    tvec{i} = (0:ntime-1)'*tscale;

    %position (condition) from which each cell came signaling data
    p{i} = s.Mats.positionIdx;
    conditions{i} = s.fixedMeta.conditions;
    
    hist_magic{i} = readmatrix(fullfile(dataDir,[dsetlabel{i} '_magic.csv']));
    hist_magic{i} = hist_magic{i}(2:end,2:end)';
    
    hist_dmagic_scramble{i} = readmatrix(fullfile(dataDir,[dsetlabel{i} '_scramble_dmagic.csv']));
    hist_dmagic_scramble{i} = hist_dmagic_scramble{i}(2:end,2:end)';
    fm_dmagic_scramble{i} = readmatrix(fullfile(dataDir,[dsetlabel{i} '_fm_scramble_dmagic.csv']));
    fm_dmagic_scramble{i} = fm_dmagic_scramble{i}(2:end,2:end)';
    
    hist_dmagic_fate{i} = readmatrix(fullfile(dataDir,[dsetlabel{i} '_dmagic_fate.csv']));
    hist_dmagic_fate{i} = hist_dmagic_fate{i}(2:end,2:end)';
    fm_dmagic_fate{i} = readmatrix(fullfile(dataDir,[dsetlabel{i} '_fm_dmagic_fate.csv']));
    fm_dmagic_fate{i} = fm_dmagic_fate{i}(2:end,2:end)';
    
    hist_dmagic_hist{i} = readmatrix(fullfile(dataDir,[dsetlabel{i} '_dmagic_hist.csv']));
    hist_dmagic_hist{i} = hist_dmagic_hist{i}(2:end,2:end)';
    fm_dmagic_hist{i} = readmatrix(fullfile(dataDir,[dsetlabel{i} '_fm_dmagic_hist.csv']));
    fm_dmagic_hist{i} = fm_dmagic_hist{i}(2:end,2:end)';
end

% redefine condition labels for single-condition dataset
di = 2;
npos = numel(unique(p{di})); nwell = 2; ppw = npos/nwell; pp = p{di}; %one condition with two wells, 3 positions per well
for ii = 1:npos
    p{di}(pp == ii) = ceil(ii/ppw);
end
conditions{di} = {[conditions{di}{1}  '_well1'], [conditions{di}{1}  '_well2']};

% redefine channel labels for repeat stains
channelLabels{di}{4} = 'DAPI2';
channelLabels{di}{8} = 'DAPI3';
channelLabels{di}{11} = 'DAPI4';
channelLabels{di}{12} = 'TFAP2C-2';
channelLabels{di}{13} = 'GATA3-2';
disp(channelLabels{di})


%% signaling density plot over time from histories
close all

highSig = cell(1,ndsets);
lowSig = cell(1,ndsets);
for di = 1:ndsets
    highSig{di} = mean(mean(hist_raw{di}(6:66,:)));% mean(meanSig(6:66));
    lowSig{di} =  mean(mean(hist_raw{di}(1:2,:)));%mean(meanSig(1:4));
end

di = 2;
npts = 100; yl = [-0.5,2.5];
pts = linspace(yl(1),yl(2),npts);
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

    meanSig(ti) = mean(data,'omitnan');
end

figure('Position',figpos); hold on
surf(T,S,X,'LineStyle','none')
cleanSubplot(fs); view(2); axis square; %clean up axes, fonts, etc
colormap(flipud(brewermap([],"RdBu"))); set(gca,'CLim',[0 1.7]); %colormap and its slimits
xlim(tvec{di}([1,end])); ylim(yl)
xlabel('time (hr)'); ylabel('SMAD4 (N:C)')
plot3(tvec{di}, meanSig, tvec{di}*0 + max(X(:)),'LineWidth',3,'Color','k')
hold off; grid off; box on; set(gca,'layer','top');
legend({'probability', 'mean'})
savefigure(fullfile(savedir,[dsetlabel{di} '_distributionVtimeFromHist']));

%% make fate marker correlation heatmap
close all

di = 2;
disp(dsetlabel{di})
lt = load(fullfile(dataDir,[dsetlabel{di}(1:6) '_lt.mat']));
lt = lt.lt;

liveMeta = lt.liveMeta;
fixedMeta = lt.fixedMeta;

channelLabel = fixedMeta.channelLabel;
npos = length(lt.live_position);

%channels in order for heatmap; no dapi, only one channel for markers for which multiple stains were done
chans = [4 3 6 13 14 2 8];
nc = length(chans);
cls = channelLabel(chans);

m = cell(npos,1);
for ii = 1:npos
    %get expression of non-dividing SMAD4 cells
    mask = lt.fixed_position(ii).cellData.labels == 1; %mask for non-dividing cells
    mask = mask & lt.mapped_idxs{ii}' > 0; %only look at fixed cells mapped to live cells
    m{ii} = lt.fixed_position(ii).cellData.nucLevel(mask,chans); %nuclear intensity in the chosen channels in these cells
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
cleanSubplot(fs); colormap(flipud(brewermap([],"RdBu")))
set(gca, 'YTick', 1:length(idxs), 'YTickLabel', cls(idxs));
set(gca, 'XTick', 1:length(idxs), 'XTickLabel', cls(idxs));
xtickangle(gca,45)
colorbar('Ticks',[-1 0 1]);
saveas(gcf,fullfile(savedir,'pairwise_correlations.png'))
    
%% fit a 2D GMM for visualization
%can't visualize a high-dimensional Gaussian mixture model so show a 2
%dimensional fit to 2 representative channels (one amnion, one pluri)
close all

k = 2; %number of components for the GMM
suffix = sprintf('_k%d',k);

%channels for which to visualize the gaussian fit (6 = isl1, 1 = nanog)
chan = [6 1];
xyz = X(:, chan);
cl = channelLabels{di}(chan);

options = statset('Display','final');
gm = fitgmdist(xyz,k,'Options',options);

figure('Position',figpos); hold on
scatter(xyz(:,1),xyz(:,2),10,'k','filled')
gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(gm,[x0,y0]),x,y);
XL = [-2 2]; YL = [-2 2]; %axis limits
fcontour(gmPDF,[XL, YL],'LineWidth',2); 
cleanSubplot(fs);% axis square
xlabel(['log( ' cls{chan(1)} ' )']); ylabel(['log( ' cls{chan(2)} ' )'])
xlim(XL); ylim(YL)
saveas(gcf,fullfile(savedir,['2d_gmm_fit',cell2mat(strcat('_',cls(chan))),suffix,'.png']))


%% plot a subset of example signaling histories
%plot 50 signaling histories with and without denoising from magic
di = 2; %dataset
disp(dsetlabel{di})

nhist = 50; %number of example histories
ridx = 1:nhist; %indices of example histories -> just the first 50 here

close all
for filtertype = 1:2
    figure('Position',figpos);hold on
    if filtertype==1
        histsample = hist_raw{di}(:,ridx);
        suffix = '';
    elseif filtertype==2
        histsample = hist_magic{di}(:,ridx);
        suffix = '_magic';
    end
    %normalize
    histsample = (histsample - lowSig{di})/(highSig{di}-lowSig{di});
    
    plot(tvec{di},histsample,'LineWidth',1) %plot individual signaling histories together
    plot(tvec{di}, meanSig, 'LineWidth',3,'Color','k') %plot mean signaling over histories
    hold off
    xlim([tvec{i}(1), tvec{di}(end)]); ylim([-0.5 2.5]);
    text(2,-0.25,['N=' num2str(nhist) '/' num2str(size(hist_raw{di},2))],'FontSize',fs)
    cleanSubplot(fs,2);%  axis square;
    xlabel('time (hr)');  ylabel('SMAD4 (N:C)');
    savefigure(fullfile(savedir,[dsetlabel{di}(1:6) '_' num2str(nhist) 'randomhistories' suffix]));
end


%% example signaling history with sigmoid fit
%show the sigmoid fit with and without denoising
close all
X = hist_raw{di};
X = (X - lowSig{di})/(highSig{di}-lowSig{di});

X_mag = hist_magic{di};
X_mag = (X_mag - lowSig{di})/(highSig{di}-lowSig{di});

visualizefits = 'on';
nhists = size(X,2);

Tm = zeros(nhists,1); Us = Tm; Ls = Tm; Ds = Tm;
            
startPoint = [0.5 1 26 1]; %starting values from which to optimize the fit
colors = lines(2);
figure('Position',figpos)
for ii = 10
    clf
    v = X(:,ii);
    t = tvec{di}(4:end); %time vector excluding time before BMP treatment
    %sigmoid function definition:
    % x(1) = high; x(2) = high - low; x(3) = duration; x(4) = sharpness of sigmoid
    fun = @(x) x(1) - x(2)./(1 + exp(-(t-x(3))/x(4)));
    
    fundiff = @(x) fun(x) - X(4:end,ii);
    fundiff_mag = @(x) fun(x) - X_mag(4:end,ii);
    x0 = startPoint;
    lb = [0 0 0 0];
    ub = [2 2 42 42];
    fitparam = lsqnonlin(fundiff, x0, lb, ub);
    fitparam_mag = lsqnonlin(fundiff_mag, x0, lb, ub);
    Us(ii) = fitparam(1); Ls(ii) = fitparam(1)- fitparam(2); Tm(ii) = fitparam(3); 
    
    hold on
    plot(tvec{di}, X(:,ii), 'LineWidth',1,'Color',colors(1,:)) %raw history
    plot(tvec{di}, X_mag(:,ii), 'LineWidth',1,'Color',colors(2,:)) %denoised history
    plot(t, fun(fitparam), 'LineWidth',3,'Color',colors(1,:)) %sigmoid fit of raw history
    plot(t, fun(fitparam_mag), 'LineWidth',3,'Color',colors(2,:)) %sigmoid fit of denoised history
    legendstr = {'data','data diffused','fit','fit diffused'};
    legend(legendstr,'FontSize',20,'Position',[0.5577 0.6687 0.3875 0.2479])

    hold off
    xlim([min(tvec{di}),max(tvec{di})])
    ylim([-0.5,2.5])
    xlabel('time (hr)'); ylabel('SMAD4 (N:C)');
    cleanSubplot(fs,2)
    savefigure(fullfile(savedir,[dsetlabel{di}(1:6) '_fit' num2str(ii)]));
    pause(0.7)
end

%% single-variable cell fate
close all

chan = [1,3]; %channels to use for the log expression ratio; 1 = ISL1, 2 = SOX2, 3 = NANOG
rthresh = zeros(1,ndsets);
ratio = cell(1,ndsets);

datasetindices = 1:ndsets;

for di = datasetindices% 4 5]

    disp(dsetlabel{di})
    if length(chan) == 2
        fm1 = fm{di}(chan(1),:);
        fm2 = fm{di}(chan(2),:);
        ratio{di} = log(fm1./fm2);
        
        rlabel =...
            strcat("log( ",channelLabels{di}{chan(1)}," / ",channelLabels{di}{chan(2)}," )");
    else
        ratio{di} = log(fm{di}(chan,:));
        rlabel = strcat("log( ",channelLabels{di}{chan}," )");
    end
    
    figure('Position',figpos)
    [f,xi] = ksdensity(ratio{di});
    dx = xi(3)-xi(2);
    hold on
%     counts = histc(ratio{di}, xi);
    counts = histcounts(ratio{di}, linspace(xi(1)-dx/2,xi(end)+dx/2,length(xi)+1));
    bar(xi, counts/sum(counts.*dx));
    plot(xi,f,'LineWidth',4)
    hold off
    xlim([-3 3]);
    xlabel(rlabel); 
    ylabel('density')
    cleanSubplot(fs,2); axis square;
    box off
    xline(rthresh(di),'LineWidth',1.5,'LineStyle','--');
    savefigure(fullfile(savedir,[dsetlabel{di} '_fateseparation.png']))
    close;
end

%% signaling heatmap sorted by cell fate

idx = cell(1,ndsets);
datasetindices = 1:ndsets;
S = cell(1,ndsets);
S_raw = cell(1,ndsets);
S_magic = cell(1,ndsets);

for filtertype = 1%:4
    for di = datasetindices
        idx{di} = p{di}>0;
        fatem = fm{di}; %expression labels
        lfm = log(fatem'); %log transformed cell fate data
        
        Xp = hist_raw{di};
        [~, S_raw{di}, ~, ~, ~, ~] = pca(Xp');
        S{di} = S_raw{di};
        
        disp(dsetlabel{di})
        smin = 0; smax = 2;

        Xp = (Xp - lowSig{di})/(highSig{di}-lowSig{di});
        ratiop = lfm(idx{di},chan(1)) - lfm(idx{di},chan(2)); %log ISL1 to NANOG expression ratio
        [B,I] = sort(ratiop); %find fate variable sorting order
        
        figure('Position',figpos)
        sortedX = Xp(:,I)'; %sort signaling histories based on the order from cell fate
        sortedX = medfilt2(sortedX,[2 2]); %median filter with a small kernal
        imagesc(tvec{di}, B ,sortedX)
        colormap(flipud(brewermap([],"RdBu")));
        cleanSubplot(fs,2); pbaspect([0.9 1 1])
        h = colorbar('Ticks',[0 2]); set(get(h,'label'),'string','SMAD4 (N:C)');
        xlabel('time (hr)'); ylabel(rlabel);
        caxis([0,2]) %color limits
        set(gca,'Ydir', 'normal')
        savefigure(fullfile(savedir,[dsetlabel{di} '_signalHeatmapByFate.png']))
    end
end


%% average signaling within fates, single condition dataset
close all

di = 2; %single condition dataset
figure('Position',figpos); hold on
X = hist_raw{di}; %raw signaling data
X = (X - lowSig{di})/(highSig{di}-lowSig{di}); %normalize between high and low values
r = ratio{di}; %cell fate variable

colors = lines(2);

marg = 0.5; %margin between negative and positive cells
plot(tvec{di}, mean(X,2),'LineWidth',lw,'Color','k')
plot(tvec{di}, mean(X(:,r>0+marg),2),'LineWidth',lw,'Color',colors(1,:))
plot(tvec{di}, mean(X(:,r<0-marg),2),'LineWidth',lw,'Color',colors(2,:))
hold off
cleanSubplot(fs,2)
xlim([min(tvec{di}),max(tvec{di})]); axis square
xlabel('time (hr)'); ylabel('SMAD4 (N:C)')
legendstr = {'combined','amnion','pluripotent'};
legend(legendstr,'Location','southwest'); 
savefigure(fullfile(savedir,[dsetlabel{di} '_avgSignalingByFate.png']))

%% average signaling within fates by position overlayed
close all

di = 1;
ps = unique(p{di}); %unique condition indices
npos = numel(ps); %number of conditions/positions
X = hist_raw{di}; %signaling data

leg = {'amnion','pluripotent'}; %legend text
%colors = 0.8*[1 0 0; 0 1 1]; %colors for the plots
colors = lines(4);
figure('Position',figpos); hold on

marg = 0.5;
for pidx = 1:npos
	lfm = log(fm{di}(:,p{di}==pidx)'); %log transformed cell fate markers
    xpos = X(:,p{di}==pidx); %signaling within the 
    
    % normalize signaling from high to low
    lowSigp = mean(mean(xpos(1:2,:)));
    highSigp = mean(mean(xpos(6:60,:)));
    xpos = (xpos - lowSigp)/(highSigp-lowSigp);
    legid = true(1,2);
    
    rpos = lfm(:,1)-lfm(:,3); %log expression ratio
    
    if sum(rpos >= rthresh(di) + marg) > 10 %don't plot average signaling in a fate with less than 10 cells
        plot(tvec{di},mean(xpos(:,rpos >= rthresh(di) + marg),2),'-','LineWidth',lw,...
            'Color',colors(1,:));
    else
        legid(1) = false; %don't include a legend entry if not displaying average signaling for this condition
    end

    if sum(rpos < rthresh(di) - marg) > 10
        plot(tvec{di},mean(xpos(:,rpos < rthresh(di) - marg),2),'LineWidth',lw,...
            'Color',colors(2,:));
        disp(['neg mean: ' num2str(mean(mean(xpos(:,rpos < rthresh(di) - marg),2)))])
    else
        legid(2) = false;
    end
    
    legend(leg,'Location','southwest')
    cleanSubplot(fs); xlim([min(tvec{di}),max(tvec{di})]); axis square
    xlabel('time (hr)'); ylabel('SMAD4 (N:C)')
end
ylim([-0.1 , 1.8])

savefigure(fullfile(savedir,[dsetlabel{di} '_avgSignalingByFateByPosition.png']))


%% separate along principle components
%for previous results, either use 210801 dataset or only position/condition
%1 (30k/BMP50) for 210827 dataset

%datasetindices = 1;

S = cell(1,ndsets);
explained = cell(1,ndsets);

S_dmagic_fate = cell(1,ndsets);

S_dmagic_hist = cell(1,ndsets);

for di = 2
    
    %find principal component scores for normalized signaling histories
    Xp = hist_raw{di}(:,idx{di});
    Xp = (Xp - lowSig{di})/(highSig{di}-lowSig{di});
    [~, S{di}, ~, ~, explained{di}, ~] = pca(Xp');
    
    % separation of signaling histories along the first 3 principle components
    for ii = 1:3
        figure('Position',figpos); hold on
        s = S{di}(:,ii);
        mu = mean(s);
        sigma = std(s);
        plot(tvec{di},mean(Xp(:,s < mu - sigma),2),'LineWidth',lw)
        plot(tvec{di},mean(Xp(:,s > mu + sigma),2),'LineWidth',lw)
        plot(tvec{di},mean(Xp,2),'LineWidth',lw, 'Color','k')
        hold off
        lstr = sprintf('x_%d',ii);
        legend([lstr,'<-\sigma'],[lstr,'>\sigma'],'mean')
        %title(sprintf('PC %d (%d%%)',ii, round(explained{di}(ii))))
        %title(sprintf('principal component %d',ii))
        cleanSubplot(fs);% axis square %cleanSubplot(fs);
        xlim([0,42]); 
        ylim([-0.5 2.5])
        xlabel('time (hr)')
        if ii == 1
            ylabel('SMAD4 (N:C)')
        end
        savefigure(fullfile(savedir,[dsetlabel{di} '_signalingAlongPC' num2str(ii),'.png']))
%         close;
    end
    
end


%% determine duration, high level, low level, integral, duration over threshold

fs = 26;

%datasetindices = 4;%2;%[1 2 4 5];
%[1 2];
close all

nftypes = 5;

allstats = cell(nftypes,ndsets);
ratio_raw = cell(1,ndsets);
ratio_dmagic_hist = cell(1,ndsets);
ratio_dmagic_fate = cell(1,ndsets);
ratio_dmagic_scramble = cell(1,ndsets);
confM = {nftypes,ndsets};
S = cell(1,ndsets);
maxaccs = cell(nftypes,ndsets);

datasetindices = 1:ndsets;

for filtertype = 1:nftypes%[1 4]%[1 4]%1:4%[1 2 3 4 5]%1:3%:3

    stats = cell(1,ndsets);
    useValidated = false; %don't have validated histories included for now

    for di = datasetindices
        %compare results with different denoising schemes, including
        %without denoising and with scrambled densoising
        if filtertype == 1
            suffix = '';
            X = hist_raw{di};
            fatem = fm{di};
            lfm = log(fatem');%./mean(fatem'));
            ratio = lfm(:,1) - lfm(:,3);
            ratio_raw{di} = ratio;
            
            Xp = X(:,idx{di});
            [~, S_raw{di}, ~, ~, ~, ~] = pca(Xp');
            S{di} = S_raw{di};
            
        elseif filtertype == 2
            suffix = '_magic';
            X = hist_magic{di};
            fatem = fm{di};
            lfm = log(fatem');%./mean(fatem'));
            ratio = lfm(:,1) - lfm(:,3);
            
            Xp = X(:,idx{di});
            [~, S_magic{di}, ~, ~, ~, ~] = pca(Xp');
            S{di} = S_magic{di};
            
        elseif filtertype == 3
            suffix = '_dmagic_hist';
            X = hist_dmagic_hist{di};
            lfm = fm_dmagic_hist{di}';
            ratio = lfm(:,1) - lfm(:,3);
            ratio_dmagic_hist{di} = ratio;
            
            Xp = X(:,idx{di});
            [~, S_dmagic_hist{di}, ~, ~, ~, ~] = pca(Xp');
            S{di} = S_dmagic_hist{di};
            
        elseif filtertype == 4
            suffix = '_dmagic_fate';
            X = hist_dmagic_fate{di};
            lfm = fm_dmagic_fate{di}';
            ratio = lfm(:,1) - lfm(:,3);
            ratio_dmagic_fate{di} = ratio;
            
            Xp = X(:,idx{di});
            [~, S_dmagic_fate{di}, ~, ~, ~, ~] = pca(Xp');
            S{di} = S_dmagic_fate{di};
            
        elseif filtertype == 5
            suffix = '_scramble';
            X = hist_dmagic_scramble{di};
            lfm = fm_dmagic_scramble{di}';
            ratio = lfm(:,1) - lfm(:,3);
            ratio_dmagic_scramble{di} = ratio;
            
            [~, S{di}, ~, ~, ~, ~] = pca(X');
        end

        % normalize signaling
        X = (X - lowSig{di})/(highSig{di}-lowSig{di});
      
        % do sigmoid fitting for each cell's signaling history        
        nhists = size(X,2);
        startPoint = [0.5 1 26 1];
        Tm = zeros(nhists,1); Us = Tm; Ls = Tm; Ds = Tm; fitint = Tm;
        tic
        for ii = 1:nhists
            t = tvec{di}(4:end); % exclude time before BMP
            %sigmoid function definition
            fun = @(x) x(1) - x(2)./(1 + exp(-(t-x(3))/x(4)));
            fundiff = @(x) fun(x) - X(4:end,ii);
            x0 = startPoint;
            lb = [0 0 0 0];
            ub = [2 2 42 42];
            fitparam = lsqnonlin(fundiff, x0, lb, ub);
            %[x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(___)
            Us(ii) = fitparam(1); Ls(ii) = fitparam(1)- fitparam(2); Tm(ii) = fitparam(3); 
            % if there is no difference between high/low then duration
            % is not defined, call it NaN
            if Us(ii) - Ls(ii) < 0.05 
                Tm(ii) = NaN;
            end
            fitint(ii) = mean(fun(fitparam)); %fint integral; this should be the same as the integral of the raw data
        end
        toc

        %collect signaling features in a single matrix of size (ncells x nfeatures)
        tmax = max(tvec{di});
        rawint = mean(X(tvec{di} < tmax,:),1)'; % is nearly identical to fitint
        stats{di} = [Tm Us Ls rawint];
        stitles = {'duration (hr)','high level','low level','integral'};
        slabels = {'duration','highlevel','lowlevel','integral'};
        nstats = size(stats{di},2);
        
        %----------------------------------------------------------------
        
        allstats{filtertype, di} = stats{di};
        
        goodstatsidx = ~any(isnan(stats{di}(idx{di},:)),2);
        
        ii = 4;
        figure('Position',figpos)
        hold on

        % also scatter raw
        lfmraw = log(fm{di}');
        rraw = lfmraw(:,1)-lfmraw(:,3);
        straw = allstats{1,di}(:,ii);
        scatter(straw, rraw,'o','MarkerEdgeColor', lines(1))

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
        cleanSubplot(fs); 
        axis square;
        savefigure(fullfile(savedir,[dsetlabel{di} '_featureVsFate_' slabels{ii} '_' suffix 'corrfate.png']))
    
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
        vlabel = " all histories";
        
        maxMIs = NaN(nstats,2);
        maxaccs{filtertype,di} = NaN(nstats,2);

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
%             confM{filtertype, sidx};
            fprintf('max accuracy = %.2g\n',accmax)
            fprintf('class balanced = %.2g\n',accbmax)
            maxaccs{di}(sidx,:) = [accbmax,vals(I2)];
        %     maxaccs(sidx,:) = [accmax,vals(I1)];
        end
    
    % DISTRIBUTIONS
    %--------------------------------------
        
        ms = 40; % marker size for data scater
        
        

        %--------------------------------------
        % HEATMAP stats vs PC
        %--------------------------------------
        figure('Position',figurePosition([640 560]))
        shortertitles = {'duration','high','low','integral'};
        if filtertype == 1
            varidx = 1:3;
        else
            varidx = 1:4;
        end
        st = stats{di}(idx{di},varidx);
%         goodstatsidx = ~any(isnan(stats{di}(idx{di},:)),2);
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
        pause(1)
        close;
    end
end

%% classification via thresholding of signaling statistics

filtertypes = {'raw','magic','dmagic_hist','dmagic_fate'};

di = 1;

for filtertype = 1
    
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

    maxMIs = NaN(nstats,2);
    maxaccs{filtertype,di} = NaN(nstats,2);

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
        maxaccs{filtertype,di}(sidx,:) = [accbmax,vals(I2)];
    end
    disp(maxaccs{filtertype,di})
end

conditional = true;
separatedists = true;
nstats = 4;

di = 1;

for sidx = 1:nstats
    figure('Position',figpos); hold on
    np = 100;
%             ydata = ratio{di}(:); 
%             xdata = stats{di}(:,sidx);

    filtertype = 1;
    fatem = fm{di};
    lfm = log(fatem');
    
    ydata = lfm(:,chan(1))-lfm(:,chan(2)); 
    xdata = allstats{filtertype, di}(:,sidx);
    
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

%     truepos = sum(sum(F.*(xx>maxaccs{filtertype,di}(sidx,2)).*(yy>ft)));
%     falsepos = sum(sum(F.*(xx>maxaccs{filtertype,di}(sidx,2)).*(yy<ft)));
%     trueneg = sum(sum(F.*(xx<maxaccs{filtertype,di}(sidx,2)).*(yy<ft)));
%     falseneg = sum(sum(F.*(xx<maxaccs{filtertype,di}(sidx,2)).*(yy>ft)));
%     M = [falseneg, truepos; trueneg, falsepos];

    disp('----------')
    truepos = sum(sum(F2.*(xx>maxaccs{filtertype,di}(sidx,2)).*(yy>ft)));
    falsepos = sum(sum(F2.*(xx>maxaccs{filtertype,di}(sidx,2)).*(yy<ft)));
    trueneg = sum(sum(F2.*(xx<maxaccs{filtertype,di}(sidx,2)).*(yy<ft)));
    falseneg = sum(sum(F2.*(xx<maxaccs{filtertype,di}(sidx,2)).*(yy>ft)));
    M2 = [falseneg, truepos; trueneg, falsepos];
    
    if conditional && true
        disttype = 'conditional';
        F = F2;
%         M = M2;
    else
        disttype = 'joint';
    end

    % override with raw data confusion matrix
    M = confM{filtertype, sidx};
    M = M(:,[2 1]);

    surf(xx,yy,F,'LineStyle','none')
    hold on 
    scatter3(xdata, ydata, xdata*0+ max(F(:)),60,'.','MarkerEdgeColor','#777777');

    % plot raw and filtered on top
    ydata = ratio_raw{di};
    xdata = allstats{1, di}(:,sidx);
    scatter3(xdata, ydata, xdata*0+ max(F(:)),30,'o','MarkerEdgeColor','#777777');


    cleanSubplot(22); view(2); axis square; colormap(flipud(brewermap([],"RdBu")))
    %xlims = {[12 26],[0.85 1.05],[0.65 0.75],[0.7 0.9],[15 35]};
    %xlim(xlims{sidx}); ylim([-3 3])
    xlim(xl); ylim(yl)
    if (sidx == 1 && ~separatedists) || (sidx == 3 && separatedists)
        ylabel(rlabel)
    end
    set(gca, 'YTick', [-1 0 1]);
    xlabel(stitles{sidx})
    
    lw = 3;
    xline(t,'--','LineWidth',lw,'Color','k');
    yline(ft,'--','LineWidth',lw);
    
    set(gca,'layer','top');
    grid off
    box on
    
%     savefigure(fullfile(savedir,[dsetlabel{di} '_signalFeaturesDistributions_' disttype '_' slabels{sidx} suffix]))

    zc = max(F(:))+0.01;
    zlim([0 zc]);

    tfs = 22;
    cb ='k';
    text(xl(1) + (t-xl(1))/2, yl(2)/2, zc, [num2str(round(M(1,1)*100)) '%'], 'Color',cb,'FontSize',tfs,'HorizontalAlignment','center','FontWeight','bold');%,'BackgroundColor','white')
    text(t + (xl(2)-t)/2, yl(2)/2, zc, [num2str(round(M(1,2)*100)) '%'], 'Color',cb,'FontSize',tfs,'HorizontalAlignment','center','FontWeight','bold')
    text(xl(1) + (t-xl(1))/2, yl(1)/2, zc, [num2str(round(M(2,1)*100)) '%'], 'Color',cb,'FontSize',tfs,'HorizontalAlignment','center','FontWeight','bold')
    text(t + (xl(2)-t)/2, yl(1)/2, zc, [num2str(round(M(2,2)*100)) '%'], 'Color',cb,'FontSize',tfs,'HorizontalAlignment','center','FontWeight','bold')

    tfs = 20;
    cnb = 'w';
    text(xl(1) + (t-xl(1))/2, yl(2)/2, zc, [num2str(round(M(1,1)*100)) '%'], 'Color',cnb,'FontSize',tfs,'HorizontalAlignment','center');
    text(t + (xl(2)-t)/2, yl(2)/2, zc, [num2str(round(M(1,2)*100)) '%'], 'Color',cnb,'FontSize',tfs,'HorizontalAlignment','center');
    text(xl(1) + (t-xl(1))/2, yl(1)/2, zc, [num2str(round(M(2,1)*100)) '%'], 'Color',cnb,'FontSize',tfs,'HorizontalAlignment','center');
    text(t + (xl(2)-t)/2, yl(1)/2, zc, [num2str(round(M(2,2)*100)) '%'], 'Color',cnb,'FontSize',tfs,'HorizontalAlignment','center');

    title(sprintf('accuracy = %.2g%%; MI = %.2g',maxaccs{di}(sidx,1),maxMIs(sidx,1)))
    hold off
    pause(1)
    
%     savefigure(fullfile(savedir,[dsetlabel{di} '_signalFeaturesDistributions_' disttype '_annot_' slabels{sidx} suffix]))
end



%% scatter feature vs fate
fs = 26;

%axim limits for each parameter
ylims = {[10 35],[0.6 1.6],[0.2 1.1],[0.3 1.35]};
scsym = {'o','.','.','.'}; %scatter plot marker types
scsiz = {50,100,100,100}; %scatter plot marker sizes

di = 1;
figure('Position',figurePosition(560,560)) 
for ii = 1:4 % which fit parameters
    corrcoeff = [0 0];
    clf
    hold on
    for filtertype = [1 4] %scatter with and without denoising
        
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
        
        y = lfm(:,1) - lfm(:,3);
        x = allstats{filtertype,di}(:,ii);
        good = ~isnan(x) & ~isnan(y);
        scatter(x,y, scsiz{filtertype}, [0.4 0.4 0.4], scsym{filtertype})
    end
    hold off
    xlabel(stitles{ii}); ylabel('log( ISL1 / NANOG )')
    xlim(ylims{ii}); ylim([-2 2])
    cleanSubplot(fs,2)
    legend({'raw','denoised'},'Location','NorthWest');
%     savefigure(fullfile(savedir,[dsetlabel{di} '_' stitles{ii} '_fate_condition_' suffix]))
end


%% scatter PC vs fitparam w/ and w/o MAGIC (fig for 220928)
di = 2;

figure('Position',figurePosition(640,560))
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
    xlabel(['principal component ' num2str(pcai)]); ylabel(stitles{ii});
    ylim(ylims{ii})
    cleanSubplot(fs,2)
    legend({'raw','denoised'},'Location','SouthEast');
    savefigure(fullfile(savedir,[dsetlabel{di} '_' stitles{ii} '_PC' num2str(pcai) '_rawVmagic']))
end


%% fitparam vs fate (SI fig for 220928)
di = 2;

%set figure size based on the number of subplots
ncols = 4+1; %number of columns
nrows = 4; %number of rows
L = 180;
L2 = L+50;
fs = 14;
figure('Position',figurePosition([L*ncols L2*nrows]))
xlims = {[10 35],[0.5 1.7],[0.2 1],[0.4 1.25]};

filtertypes = [1 3 4 5];
for ii = 1:length(filtertypes) %iterate over denoising schemes
        
    for param = 1:4 %iterate over signaling features/parameters

        
        filtertype = filtertypes(ii);
        
        if filtertype <= 2
            r = ratio_raw{di};
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
    hist(r,linspace(-2,2,50)) %#ok<HIST>
    h = findobj(gca,'Type','patch');
    h.FaceColor = lines(1);
    h.EdgeColor = 'w';
    xlabel('log( ISL1 / NANOG )');
    ylabel('cell number');
    cleanSubplot(fs,lw)
    axis square
end
savefigure(fullfile(savedir,[dsetlabel{di} '_featuresVsFateAllScatters']))


%% scatter feature vs fate colored for conditions

di = 1;
figure('Position',figurePosition(560,560)) 
for ii = 1:4 % which fit parameters
    ylims = {[10 35],[0.6 1.6],[0.2 1.1],[0.3 1.35]};
    scsym = {'o','.','.','.'};
    scsiz = {50,100,100,100};
    pcai = 1; % which principal component
    corrcoeff = [0 0];
    clf
    hold on
    for filtertype = [1 4] %scatter with and without denoising
        
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
        
        y = lfm(:,1) - lfm(:,3);
        x = allstats{filtertype,di}(:,ii);
        good = ~isnan(x) & ~isnan(y);
        colormap lines
        scatter(x,y, scsiz{filtertype}, p{di}, scsym{filtertype})
        [corrcoeff(filtertype),pval] = corr(x(good),y(good));
        if filtertype>1
            text(ylims{ii}(2), -1.5, ['corr=' num2str(corrcoeff(filtertype),2)],'FontSize',20,'HorizontalAlignment','right');
            text(ylims{ii}(2), -1.75, ['p=' num2str(pval,2)],'FontSize',20,'HorizontalAlignment','right');
        end
    end
    hold off
    xlabel(stitles{ii}); ylabel('log( ISL1 / NANOG )')
    xlim(ylims{ii}); ylim([-2 2])
    cleanSubplot(fs,2)
    legend({'raw','denoised'},'Location','NorthWest');
    savefigure(fullfile(savedir,[dsetlabel{di} '_' stitles{ii} '_fate_condition_' suffix]))
end


%% variance explained with and without magic

di = 2; %single condition dataset

%pca on raw histories
Xp = hist_raw{di}(:,idx{di});
Xp = (Xp - lowSig{di})/(highSig{di}-lowSig{di});
[~, S{di}, ~, ~, explained_raw, ~] = pca(Xp');

%pca on denoised histories
Xp = hist_magic{di}(:,idx{di});
Xp = (Xp - lowSig{di})/(highSig{di}-lowSig{di});
[~, S_magic{di}, ~, ~, explained_magic, ~] = pca(Xp');

colors = lines(2);
nPC = 10; %number of principal components to plot; don't need to show a long tail of near-zero values
figure('Position',figurePosition(560,560)) ; hold on
%plot values of variance explained with and without denoising
plot(1:nPC, explained_raw(1:nPC), '-','Color',colors(1,:),'LineWidth',3) %line plot
plot(1:nPC, explained_magic(1:nPC), '-','Color',colors(2,:),'LineWidth',3)
scatter(1:nPC, explained_raw(1:nPC), 1000, '.','MarkerEdgeColor',colors(1,:)) %mark discrete values
scatter(1:nPC, explained_magic(1:nPC), 1000, '.','MarkerEdgeColor',colors(2,:))
legend({'raw','denoised'})
xlabel('principal component'); ylabel('var explained (%)');
hold off
cleanSubplot(fs,2); axis square
savefigure(fullfile(savedir,[dsetlabel{di} '_varianceExplained.png']))










