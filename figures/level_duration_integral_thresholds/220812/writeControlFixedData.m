clear; close all; clc

%%

scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
dataDir = scriptPath;
loadDir = 'Z:\Heemskerk_Lab\Seth-06\220812_Smad4GFP_BMP_IWP2_LDN_durations';

test = load(fullfile(loadDir,'combinedExpressionData.mat'));

cond1 = [1,10];
disp(test.conditions(cond1))

cond2 = 9;
disp(test.conditions(cond2))

m1 = test.m(cond1,:);
m1 = cell2mat(m1(:));

m2 = test.m(cond2,:);
m2 = cell2mat(m2(:));

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




