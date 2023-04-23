clear; close all; clc

%% set up
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
dataDir = scriptPath;
savedir = fullfile(dataDir,'figures');
if ~exist(savedir,'dir'), mkdir(savedir); end

fs = 28;
lfs = 20;
lw = 3;
scale = 0.95;

xx = linspace(0,43,500);
yy = linspace(0,1,500);

dthresh = 25;
L0 = 0.2;
ithresh = dthresh*1 + (42 - dthresh)*L0;
lline = (ithresh - 42*L0)./xx + L0;

levels = [1 0.7 0.78];
durations = [27 42 27];

fwh = [560 560]; %figure width and height
fposition = figurePosition(fwh);


%% example signaling history to show features

tau = 1;
tvec = linspace(-1,42,100);

L = 0.76;
L0 = 0.15;
D = 25;

v = signalingDynamics(tvec,L,D,tau,0);

figure('Position',fposition)
plot(tvec,v,'LineWidth',lw)
cleanSubplot(fs); axis square
set(gca,'Box','off')
xlim(tvec([1,end])); ylim([-0.05,1.05])

hold on
fill(tvec([1:end,end,1]),[v(1:end),0,0],[0 0.4470 0.7410],'FaceAlpha',0.5,'EdgeColor','none')
drawnow

xlabel('time (hr)'); ylabel('signaling (au)')
title({'','signaling features'})

pos = get(gca,'Position');
pos(2) = pos(2) + (1 - scale)*pos(4);
pos(4) = scale*pos(4);
set(gca,'Position',pos)
savefigure(fullfile(savedir,'exampleSignalingWithLabels'))

%% compare duration and level thresholding
textx = [31 39 20]; texty = [1 0.62 0.85]; textt = ["1","2","3"];

lthresh = 0.74;

%make example signaling histories
v = cell(length(durations));
for ii = 1:length(durations)
    v{ii} = signalingDynamics(tvec,levels(ii),durations(ii));
end

titles = {{'level + duration', 'thresholds'},{'','integral threshold'}};
savenames = {'durationThresholded','integralThresholded'};
savenames = fullfile(savedir,savenames);
cols = 0.8*[1 0 0; 0 1 1];
% colids = [1 1 1; 1 1 2];
colids = [2 1 2; 2 2 1];

for ii = 1:length(titles)
    figure('Position',fposition); hold on
    %plot and label signaling histories colored for fate
    for jj = 1:length(v)
        plot(tvec,v{jj},'LineWidth',lw,'Color',cols(colids(ii,jj),:))
%         text(textx(jj),texty(jj),textt(jj),'FontSize',fs,'FontWeight','bold')
    end
    %draw thresholds
    xline(dthresh,'LineStyle','--','Color','k','LineWidth',lw);
    yline(lthresh,'LineStyle','--','Color','k','LineWidth',lw);
    
    cleanSubplot(fs); axis square
    xlim(tvec([1,end])); ylim([-0.05,1.05])
    xlabel('time (hr)'); ylabel('signaling (au)')
    title(titles{ii})
    
    if ii == 1
        legend('amnion','pluripotent','Location','southwest','FontSize',lfs)
    end
    
    pos = get(gca,'Position');
    pos(2) = pos(2) + (1 - scale)*pos(4);
    pos(4) = scale*pos(4);
    set(gca,'Position',pos)
    
    savefigure(savenames{ii})
end


%% functions

function y = signalingDynamics(tvec,L,D,tau,L0)

if ~exist('tau','var')
    tau = 1;
end
if ~exist('L0','var')
    L0 = 0;
end

zeros(size(tvec));
y(tvec >= 0) = L*(1 - exp(-tau*tvec(tvec >= 0)));
y(tvec >= D) = (L-L0)*(exp(-tau*(tvec(tvec >= D) - D))+L0);

end

