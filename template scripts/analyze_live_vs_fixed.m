clear all; close all;
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);

% raw data location, modify if not the same as the location of this script
dataDir = scriptPath; 
cd(dataDir);

liveDir = fullfile(dataDir, '200114_ERKgeminin_live_inhibitorGF_thenFix_2');
fixedDir = fullfile(dataDir, '200116_ERKgeminin_live_inhibitorGF_thenFix_2_pERK_Z');

liveMeta = load(fullfile(liveDir,'meta.mat'));
liveMeta = liveMeta.meta;

fixedMeta = load(fullfile(fixedDir,'meta.mat'));
fixedMeta = fixedMeta.meta;

% %%
% positions_live = load(fullfile(liveDir,'positions.mat'));
% postions_fixed = load(fullfile(fixedDir,'positions.mat'));

pERKchannel = 1;
ERKKTRchannel = 1;
pERKchannelIndex = 2;

%% map live and fixed

lt = LineageTrace(liveDir, fixedDir);

maxDist = 10;
lt.mapPoints(maxDist);

%% check mapping between live and fixed

pi = 3;
fixedIm = lt.fixed_position(pi).loadImage(fixedDir, lt.fixed_position(pi).nucChannel, 1);
liveIm = lt.live_position(pi).loadImage(liveDir, lt.live_position(pi).nucChannel, lt.live_position(pi).nTime);

% first image is green
imshowpair(imadjust(max(fixedIm,[],3)), imadjust(liveIm));

%liveImERK = lt.live_position(pi).loadImage(liveDir, ERKKTRchannel, lt.live_position(pi).nTime);
%imshowpair(imadjust(liveImERK), imadjust(liveIm));

fixedXY = lt.fixed_position(pi).cellData(1).XY;
liveXY = lt.live_position(pi).cellData(end).XY;

mappedSourceIdx = find(lt.mapped_idxs{pi} > 0);
mappedTargetIdx = lt.mapped_idxs{pi}(lt.mapped_idxs{pi} > 0);
dXY = liveXY(mappedTargetIdx,:) - fixedXY(mappedSourceIdx,:);

hold on
scatter(fixedXY(:,1),fixedXY(:,2),10,'filled','MarkerFaceColor','blue');
scatter(liveXY(:,1),liveXY(:,2),10,'filled','MarkerFaceColor','y');
quiver(fixedXY(mappedSourceIdx,1),fixedXY(mappedSourceIdx,2),dXY(:,1), dXY(:,2),0,'LineWidth',2);

meandXY = median(dXY);
quiver(512,512,meandXY(:,1),meandXY(:,2),0,'LineWidth',2)

hold off

saveas(gcf,sprintf('fixed_live_mapping_p%.4d.png', pi));

%% look at mean trends in live v fixed

nWells = liveMeta.nWells;
posPerCondition = liveMeta.posPerCondition;
conditions = {'tesr', 'F/C', 'Mi0.5/C', 'Mi5/C','Mi5','Mi0.5', 'F','tesr'};

clf 
hold on
for wellnr = 1:nWells

    conditionPositions = posPerCondition*(wellnr-1)+1:posPerCondition*wellnr;
    
    for pi = conditionPositions
        
        liveC = lt.live_position(pi).cellData(end).cytLevel;
        liveN = lt.live_position(pi).cellData(end).nucLevel;
        liveBG = lt.live_position(pi).cellData(end).background;
        lt.live_position(pi).cellData(end).CNratio = (liveC - liveBG)./(liveN - liveBG);
        
        R = nanmean(lt.live_position(pi).cellData(end).CNratio);
        eR = nanstd(lt.live_position(pi).cellData(end).CNratio);
        errorbar(wellnr, R, eR, '.', 'MarkerSize',10);
    end
end
hold off

xticklabels(conditions)
ylim([0 2]);
saveas(gcf, fullfile(dataDir, 'live_final.png'));

clf 
hold on
for wellnr = 1:nWells

    conditionPositions = posPerCondition*(wellnr-1)+1:posPerCondition*wellnr;
    
    for pi = conditionPositions
        fixedC = lt.fixed_position(pi).cellData(1).cytLevel(:, pERKchannelIndex);
        fixedCmean = mean(fixedC);
        fixedCstd = std(fixedC);
        errorbar(wellnr, fixedCmean, fixedCstd, '.', 'MarkerSize', 10);
    end
end
hold off

xticklabels(conditions)
ylim([2000 8000]);
saveas(gcf, fullfile(dataDir, 'fixed_final.png'));

%% scatter fixed pERK vs live ERK-KTR

for pi = 3

    figure, 
    mappedSourceIdx = find(lt.mapped_idxs{pi} > 0);
    mappedTargetIdx = lt.mapped_idxs{pi}(lt.mapped_idxs{pi} > 0);

    ERKKTRvals = lt.live_position(pi).cellData(end).CNratio(mappedTargetIdx, ERKKTRchannel);
    ERKKTRvals = ERKKTRvals(~isnan(ERKKTRvals));
    pERKvals = lt.fixed_position(pi).cellData(1).cytLevel(mappedSourceIdx(~isnan(ERKKTRvals)), pERKchannelIndex);

    scatter(pERKvals, ERKKTRvals)
    corrcoef(pERKvals, ERKKTRvals)
    hold on
    scatter(mean(pERKvals),mean(ERKKTRvals),'filled','r')
    hold off
    
    ylim([0.5 3]);
    xlim([2000 10000]);
    ylabel('ERK-KTR (C:N)');
    xlabel('pERK (au)');
    axis square
    saveas(gcf, sprintf(fullfile(dataDir, 'fixed_live_correlation_p%.4d.png'), pi));
end

