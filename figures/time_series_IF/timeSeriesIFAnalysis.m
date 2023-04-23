clear; close all; clc

%%
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
dataDir = scriptPath;
savedir = fullfile(dataDir,'figures');

fs = 28; lw = 3; lfs = 20;

dirs = {'round1','round2','round3','round4'};

times = [0 6 12 24 36 48];
doses = [0 10 30 100 300];

load(fullfile(dataDir,'live','signalData.mat'))
%get average level in each dose between 0 and 24 hours and normalize
%between 0 and 1
SMAD4 = mean(X(5:146,:),1);
SMAD4 = (SMAD4 - min(SMAD4))/(max(SMAD4) - min(SMAD4));

nd = length(dirs); %number of directories
Y = cell(1,1,nd); %expression data
Err = cell(1,1,nd); %standard deviation for error bars
cls = cell(1,nd); %channel names
m = cell(1,nd); %metadata

for di = 1:nd
    cdata = load(fullfile(dataDir,dirs{di},'expressionMatrix.mat'));
    Y{di} = cdata.Y;
    Err{di} = cdata.Err;
    cls{di} = cdata.cls;
    m{di} = load(fullfile(dataDir,dirs{di},'meta.mat'));
    m{di} = m{di}.meta;
end

Y = cell2mat(Y);
Err = cell2mat(Err);
cls = cat(2,cls{:});

ndose = size(Y,2); %number of doses (treatment conditions)

%% plot SMAD4 level vs LDN dose
%can't exactly make the axis proportional given that negative control does
%not have a simple numerical equivalent LDN dose and positive control (no
%LDN) requires an offset to make sense on a log scale

figure('Position',figurePosition(560,560))

xs = log10(1 + [doses, 900]);% xs = 1:length(SMAD4);
xs(1) = 0.7;
% xs = [1 1.75 2.75 3.875 5.125 6.875];

plot(xs,SMAD4,'x','LineWidth',lw,'MarkerSize',15)
cleanSubplot(fs); axis square; xlim(xs([1,end]) + [-0.25 0.5])
set(gca,'Box','off')
set(gca,'XTick',xs,'XTickLabel',{'0','10','30','100','300','X'})
xlim([min(xs),max(xs)])
% xtickangle(30)

xlabel('BMPRi (nM)')
ylabel('SMAD4 level')
savefigure(fullfile(savedir,'SMAD4_level_vs_LDN.png'))

%% plot individual genes against time
close all

gidxs = [1 2 3 4 8 10 12];
legidxs = [1 4];
legstr = {'0','10','30','100','300','X'};
colors = turbo(ndose); colors = colors(end:-1:1,:);

for ii = 1:length(gidxs)
    gidx = gidxs(ii);
    cl = cls{gidx};
    figure('Position',figurePosition(560,560)); hold on
    for jj = 1:ndose
        ydata = Y(:,jj,gidx);
        err = Err(:,jj,gidx);
        mask = ~isnan(ydata);
        errorbar(times(mask),ydata(mask),err(mask),...
            'LineWidth',lw,'Color',colors(jj,:))
    end
    cleanSubplot(fs); xlabel('time (hr)'); ylabel(strcat(cl," (au)"))
    xlim([0,48])
    
    if ismember(gidx, legidxs)
        lgd = legend(legstr,'FontSize',lfs,'Position',[0.7423 0.6039 0.2518 0.3330]);
        title(lgd,'BMPRi (nM)','FontSize',lfs)
        disp('adjust legend')
        pause
    end
    savefigure(fullfile(savedir,['expressionVTime_',cl]))
end



%% plot expression with respect to signaling level for each time

tidxs = [3 3 6];
gidxs = {[2,8],[3,4],[1,12,10]};

for jj = 1:length(tidxs)
    ti = tidxs(jj);
    time = times(ti); tstr = num2str(time);
    cl = cls(gidxs{jj});
    figure('Position',figurePosition(560,560)); hold on
    for ii = 1:length(gidxs{jj})
        gidx = gidxs{jj}(ii);
        errorbar(SMAD4,Y(ti,:,gidx),Err(ti,:,gidx),...
            'LineWidth',lw)
    end
    cleanSubplot(fs); axis square
    set(gca,'Box','off')
    xlabel('SMAD4 level'); ylabel("intensity (au)")
    legend(cl,'FontSize',lfs)
    disp('adjust legend')
    pause
    c = strcat('_',cl);
    savefigure(fullfile(savedir,['expressionVSMAD4',cat(2,c{:}),'_',tstr,'hr.png']))
end









