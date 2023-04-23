clear; close all; clc

%%

scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
dataDir = scriptPath;
loadDir = 'Z:\Heemskerk_Lab\Seth-06\220629_Smad4GFP_BMP_IWP2_LDN_live\fixed_new\timeSeries';

% test = load(fullfile(loadDir,'combinedExpressionData.mat'));

load(fullfile(loadDir,'meta.mat'))
load(fullfile(loadDir,'positions.mat'))

ncond = length(meta.conditions);
npos = length(positions);
ppc = npos/ncond;

M = cell(ncond,1);
for cidx = 1:ncond
    m = cell(ppc,1);
    for cpi = 1:ppc
        pidx = (cidx - 1)*ppc + cpi;
        m{cpi} = positions(pidx).cellData.nucLevel;
    end
    M{cidx} = cell2mat(m);
end


%%
cond1 = 10;
disp(meta.conditions(cond1))

cond2 = 4;
disp(meta.conditions(cond2))

m1 = M(cond1,:);
m1 = cell2mat(m1(:));

m2 = M(cond2,:);
m2 = cell2mat(m2(:));

%%
m = {m1,m2};
conditions = {'mTeSR','LDN0,42hr'};

% mall = [m1;m2];
% ISL1 = mall(:,2); NANOG = mall(:,4);

chan = [2 4];

figure('Position',figurePosition(560,560)); hold on
mlog = cell(size(m));
ratios = cell(size(m));
for ii = 1:2
    xy = m{ii}(:,chan);
    xy(xy < 0) = 1;
    xy = log(xy + 1);
    mlog{ii} = xy;
    ratios{ii} = xy(:,1) - xy(:,2);
    
    [f,xi] = ksdensity(ratios{ii});
    plot(xi,f,'LineWidth',3)
end
cleanSubplot(24); axis square
xlabel('log(ISL1 / NANOG)')
ylabel('density')


save(fullfile(dataDir,'ControlConditionsIF.mat'),'m','conditions')




