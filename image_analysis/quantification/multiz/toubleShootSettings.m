clear; close all; clc

%% setup, load images & segmentations
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
dataDir = scriptPath;

%'stitched_p%.4d_w%.4d_t%.4d'
fbare = 'stitched_p%.4d_w%.4d_t%.4d';
mipbare = 'stitched_MIP_p%.4d_w%.4d_t%.4d';

nucChannel = 0;
fgchannel = 1;

pidx = 11;
ti = 140;

cytsegname = [sprintf(fbare,pidx-1,fgchannel,ti-1),'_Simple Segmentation.h5'];
cytsegname = fullfile(dataDir,cytsegname);
fgseg = ilastikRead(cytsegname);
m = size(fgseg,1); n = size(fgseg,2); nz = size(fgseg,3);

nucsegname = [sprintf(fbare,pidx-1,nucChannel,ti-1),'_FinalSegmentation.tif'];
nucsegname = fullfile(dataDir,nucsegname);
nucseg = false(m,n,nz);
for zi = 1:nz
    nucseg(:,:,zi) = imread(nucsegname,zi) == 1;
end

nucname = [sprintf(fbare,pidx-1,nucChannel,ti-1),'.tif'];
nucname = fullfile(dataDir,nucname);

nucim = zeros(m,n,nz,'uint16');
for zi = 1:nz
    fprintf('.');
    nucim(:,:,zi) = imread(nucname,zi);
end
nucmip = max(nucim,[],3);
nuclim = seglim(nucmip,sum(nucseg,3)>0);

fgname = [sprintf(fbare,pidx-1,fgchannel,ti-1),'.tif'];
fgname = fullfile(dataDir,fgname);
fgim = zeros(m,n,nz,'uint16');
for zi = 1:nz
    fprintf('.');
    fgim(:,:,zi) = imread(fgname,zi);
end
fgmip = max(fgim,[],3);
fglim = seglim(fgmip,sum(fgseg,3)>0);
fprintf('\n')

%%
zi = 3;
opts = struct(...
                    'dataChannels',     0:1,...
                    'cytoplasmicLevels',true,...
                    'fgChannel',        fgchannel,...
                    'cytoSize',         4,...
                    'cytoMargin',       0,...
                    'segmentationDir',  fullfile(dataDir,'MIP'),...
                    'MIPidxDir',        fullfile(dataDir,'MIP'),...
                    'nucShrinkage',     2,...
                    'bgMargin',         4,...
                    'segFG',            1);
imgs = {nucim(:,:,zi),fgim(:,:,zi)};
[debugInfo, cellData] =...
    extractDataSingleTime(imgs,nucseg(:,:,zi),fgseg(:,:,zi),opts);

fprintf('cytoSize = %d\n',opts.cytoSize)
fprintf('cytoMargin = %d\n',opts.cytoMargin)
fprintf('nucShrinkage = %d\n',opts.nucShrinkage)
fprintf('bgMargin = %d\n\n',opts.bgMargin)

fprintf('mean nuc:cyt ratio = %g\n mean nuc = %g\nmean cyt = %g\nbackground = %g\n\n',...
    mean(cellData.NCratio(:,2),'omitnan'),mean(cellData.nucLevel(:,2),'omitnan'),...
    mean(cellData.cytLevel(:,2),'omitnan'),cellData.background(2))

disp('--------------------')

%% interactively overlay entire nuclear and cytoplasmic masks
A = repmat(im2double(imadjust(nucim(:,:,zi),nuclim)),1,1,3);
B = repmat(im2double(imadjust(fgim(:,:,zi),fglim)),1,1,3);

s = 0.3;

overlays = cat(3, debugInfo.nucmaskraw, debugInfo.nucmask, debugInfo.cytmask);
nucoverlay = (1-s)*A + s*overlays;
cytoverlay = (1-s)*B + s*overlays;

f = figure('WindowState','maximized');
ax1 = subplot_tight(1,2,1);
p1 = imshow(nucoverlay);
cleanSubplot
title('\color{yellow}nucleus \color{blue}cytoplasm \color{red}removed nucleus')

ax2 = subplot_tight(1,2,2);
p2 = imshow(cytoverlay);
cleanSubplot
title(sprintf('mean nuc:cyt ratio = %g',mean(cellData.NCratio(:,2))))
linkaxes([ax1,ax2])
drawnow

criterion = true;
while criterion
    waitforbuttonpress;
    key = f.CurrentCharacter;
    if strcmp(key,'o')
        %show overlays
        set(p1,'CData',nucoverlay)
        set(p2,'CData',cytoverlay)
    elseif strcmp(key,'i')
        %show images
        set(p1,'CData',A)
        set(p2,'CData',B)
    elseif strcmp(key,'+')
        %zoom in 2X
        xl = get(gca,'xlim');
        yl = get(gca,'ylim');
        set(gca,'xlim', mean(xl) + 0.25*(xl(2) - xl(1))*[-1 1])
        set(gca,'ylim', mean(yl) + 0.25*(yl(2) - yl(1))*[-1 1])
    elseif strcmp(key,'-')
        %zoom out 2X
        xl = get(gca,'xlim');
        yl = get(gca,'ylim');
        set(gca,'xlim', mean(xl) + (xl(2) - xl(1))*[-1 1])
        set(gca,'ylim', mean(yl) + (yl(2) - yl(1))*[-1 1])
    elseif strcmp(key,'2')
        %pan down
        yl = get(gca,'ylim');
        set(gca,'ylim', yl + 0.25*(yl(2) - yl(1)))
    elseif strcmp(key,'8')
        %pan up
        yl = get(gca,'ylim');
        set(gca,'ylim', yl - 0.25*(yl(2) - yl(1)))
    elseif strcmp(key,'4')
        %pan left
        xl = get(gca,'xlim');
        set(gca,'xlim', xl - 0.25*(xl(2) - xl(1)))
    elseif strcmp(key,'6')
        %pan right
        xl = get(gca,'xlim');
        set(gca,'xlim', xl + 0.25*(xl(2) - xl(1)))
    elseif strcmp(key,'5')
        %reset field of view to default (full image)
        set(gca,'xlim', [0,size(A,1)] + [-0.5 0.5])
        set(gca,'ylim', [0,size(A,2)] + [-0.5 0.5])
    elseif strcmp(key,'e')
        %exit the while loop
        criterion = false;
    end
    
end

close(f)

%% read in and process images in parallel
treatmentTime = 4;
tmax = 15;
tscale = 10/60;
tvec = (1-treatmentTime:tmax-treatmentTime)*tscale;

opts = struct(...
                    'dataChannels',     0:1,...
                    'cytoplasmicLevels',true,...
                    'fgChannel',        fgchannel,...
                    'cytoSize',         4,...
                    'cytoMargin',       -1,...
                    'segmentationDir',  fullfile(dataDir,'MIP'),...
                    'MIPidxDir',        fullfile(dataDir,'MIP'),...
                    'nucShrinkage',     3,...
                    'bgMargin',         4,...
                    'segFG',            1);

clear cellData

nucsegname = sprintf(segbare,pidx-1);
nucsegname = fullfile(dataDir,'MIP',nucsegname);

cytsegname = [sprintf(mipbare,pidx-1,fgchannel),'_Simple Segmentation.h5'];
cytsegname = fullfile(dataDir,'MIP',cytsegname);
Fgseg = ilastikRead(cytsegname);

fname = sprintf(fbare,pidx-1);
fname = fullfile(dataDir,fname);
r = bfGetReader(fname);
nz = r.getSizeZ; m = r.getSizeY; n = r.getSizeX;

tic
cdata = cell(1,tmax);
for ti = 1:tmax
    nucseg = imread(nucsegname,ti) == 1;
    fgseg = Fgseg(:,:,ti);
    
    nucim = zeros(m,n,nz,'uint16');
    for zi = 1:nz
        nucim(:,:,zi) = bfGetPlane(r, r.getIndex(zi-1,nucChannel,ti-1)+1);
    end

    fgim = zeros(m,n,nz,'uint16');
    for zi = 1:nz
        fgim(:,:,zi) = bfGetPlane(r, r.getIndex(zi-1,fgchannel,ti-1)+1);
    end
    
    imgs = {nucim,fgim};
    [debugInfo, cellData] = extractDataSingleTime(imgs,nucseg,fgseg,opts);
    cdata{ti} = cellData;
end
fields = fieldnames(cdata{1});
cellData = struct;
for ti = 1:length(cdata)
    for fi = 1:length(fields)
        cellData(ti).(fields{fi}) = cdata{ti}.(fields{fi});
    end
end
toc

v = cellfun(@(x) median(x(:,2),'omitnan'),{cellData.NCratio}');
figure
plot(tvec,v,'LineWidth',2.5)
xlim(tvec([1,tmax]))
cleanSubplot(16)
xlabel('time (hours)')
ylabel('smad4 nuc:cyt ratio')



