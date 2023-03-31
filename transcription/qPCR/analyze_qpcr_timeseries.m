clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 
dataDir = '/Users/idse/data_tmp/qpcr/1702_qpcrS4Aseries';

filenames = {'170213 Smad4 ATP5O CER LHX1_data.txt'};

% filenames = {'170115 A New TimeSeries Mixl1 Pai1 Snai1 Skil_data.txt',...
%             '170115 A New TimeSeries Mixl1 Pai1 Snai1Skil B_data.txt'};
% 
% filenamesNorm = {'170115 A New TimeSeries ATP5 Oct4 Nanog Sox2_dataX.txt',...
%     '170115 A New TimeSeries ATP5O Oct4 Nanog Sox2 B_data.txt'};

normName = 'ATP5O';

%% read data

tolerance = 1;
data = combineQPCR(dataDir, filenames, tolerance);

Ntargets = data.Ntargets;
Nsamples = data.Nsamples;
targets = data.targets;
CTmean = data.CTmean;
CTstd = data.CTstd;

% normalization
refIdx = find(strcmp(targets,normName));
if isempty(refIdx)
    disp('using normalization from other file');
	refData =  combineQPCR(dataDir, filenamesNorm);
    refIdx = find(strcmp(refData.targets,normName));
    CTref = refData.CTmean(:,refIdx);
    refIdx = 0; % to make sure this idx isn't excluded from plot later
else
    disp('using normalization from same file');
    CTref = CTmean(:,refIdx);
end

barefname = [targets{:}];

%% normalize

CTnorm = zeros([Nsamples Ntargets]);
for targeti = 1:Ntargets
	CTnorm(:,targeti) = -(CTmean(:,targeti) - CTref);
    CTnorm(:,targeti) = CTnorm(:,targeti) - CTnorm(1,targeti);
end

%% raw value

Atimes = [0 1 2 4 6 9 12 24];
Btimes = [0 24];
SBtimes = 12;

Aidx = [1 3:5 7:10];
Bidx = [1 2];
SBidx = 6;

clf 
bar(Atimes, CTmean(Aidx,:))
legend(targets,'Location','SouthEast')
ylim([20 35]);
saveas(gcf, fullfile(dataDir,['CTraw_A_' barefname]));
saveas(gcf, fullfile(dataDir,['CTraw_A_' barefname '.png']));

figure,
bar(Btimes, CTmean(Bidx,:))
legend(targets,'Location','SouthEast')
ylim([20 35]);
saveas(gcf, fullfile(dataDir,['CTraw_B_' barefname]));
saveas(gcf, fullfile(dataDir,['CTraw_B_' barefname '.png']));

% figure,
% bar(SBtimes, CTmean(SBidx,:))
% legend(targets,'Location','SouthEast')
% ylim([20 35]);

%% normalized

% bar(times, CTnorm(:,2:end))
% legend(D.targets(2:end))
figure,

targeti = setdiff(1:Ntargets, refIdx);

clf
color = lines(numel(targeti));
hold on
for i = 1:numel(targeti)
    plot(Atimes, CTnorm(Aidx,targeti(i)),...
       '-x','Color',color(i,:),'LineWidth',2);
end

legend(targets(targeti),'location','SouthEast')

for i = 1:numel(targeti)
    plot(Btimes, CTnorm(Bidx,targeti(i)),...
       '--x','Color',color(i,:),'LineWidth',2);
    plot(SBtimes, CTnorm(SBidx,targeti(i)),...
       'o','Color',color(i,:),'LineWidth',2);

    errorbar(Atimes, CTnorm(Aidx,targeti(i)), CTstd(Aidx,targeti(i)),...
        '-x','Color',color(i,:),'LineWidth',2);
    errorbar(Btimes, CTnorm(Bidx,targeti(i)), CTstd(Bidx,targeti(i)),...
        '-x','Color',color(i,:),'LineWidth',2);
    errorbar(SBtimes, CTnorm(SBidx,targeti(i)), CTstd(SBidx,targeti(i)),...
        '-x','Color',color(i,:),'LineWidth',2);
end

hold off
xlim([0 24]);
%ylim([-0.5 5.5]);


fs = 15;
xlabel('time (hrs)', 'FontSize',fs, 'FontWeight','Bold');
ylabel('C_T - C_{T,ATP5O}', 'FontSize',fs, 'FontWeight','Bold');
set(gcf,'color','w');
set(gca, 'LineWidth', 2);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')

saveas(gcf, fullfile(dataDir, ['CTnorm_' barefname]));
saveas(gcf, fullfile(dataDir, ['CTnorm_' barefname '.png']));

%% normalized exponential

figure,

targeti = setdiff(1:Ntargets, refIdx);

clf
color = lines(numel(targeti));
hold on
for i = 1:numel(targeti)
    plot(Atimes, 2.^CTnorm(Aidx,targeti(i)),...
       '-x','Color',color(i,:),'LineWidth',2);
end

legend(targets(targeti),'location','NorthEast')

for i = 1:numel(targeti)
    plot(Btimes, 2.^CTnorm(Bidx,targeti(i)),...
       '--x','Color',color(i,:),'LineWidth',2);
    plot(SBtimes, 2.^CTnorm(SBidx,targeti(i)),...
       'o','Color',color(i,:),'LineWidth',2);
end

hold off
xlim([0 24]);
%ylim([-0.5 5.5]);


fs = 15;
xlabel('time (hrs)', 'FontSize',fs, 'FontWeight','Bold');
ylabel('2.\^(C_T - C_{T,ATP5O})', 'FontSize',fs, 'FontWeight','Bold');
set(gcf,'color','w');
set(gca, 'LineWidth', 2);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')

saveas(gcf, fullfile(dataDir, ['expnorm_' barefname]));
saveas(gcf, fullfile(dataDir, ['expnorm_' barefname '.png']));

    