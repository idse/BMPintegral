clear; close all; clc

%% setup
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
dataDir = scriptPath;
mipDir = fullfile(dataDir,'MIP');
bare = 'stitched_MIP_p%.4d_w%.4d_t0000.tif';

conditions = {'160k,BMP+IWP2','160k,BMP','240k,BMP','240k,BMP+IWP2'};

ncond = length(conditions);
ppc = 5;
npos = ncond*ppc;
channelLabel = {'DAPI','ISL1','pSmad1','TBXT'};
nc = length(channelLabel);

if ppc > ncond
    gridsize = [ncond ppc];
    layout = [kron((1:ncond)',ones(ppc,1)), repmat((1:ppc)',ncond,1)];
else
    gridsize = [ppc ncond];
    layout = [repmat((1:ppc)',ncond,1), kron((1:ncond)',ones(ppc,1))];
end

% gridsize = [4 4];
% layout = NaN(npos,2);
% layout(1:3:npos,:) = [kron((1:4)',ones(4,1)), repmat((1:4)',4,1)];

%% load images
imgs = cell(gridsize(1),gridsize(2),nc);
for pidx = 1:npos
    ii = layout(pidx,1); jj = layout(pidx,2);
    if ~isnan(ii)
        for ci = 1:nc
            fprintf('.')
            fname = fullfile(mipDir,sprintf(bare,pidx-1,ci-1));
            im = imread(fname);
            imgs{ii,jj,ci} = im;
        end
    end
    fprintf('\n')
end
I = find(cellfun(@(x) ~isempty(x), imgs(:,:,1)),1);
[row,col] = ind2sub(size(imgs,[1,2]),I);
mn = size(imgs{row,col,1}); m = mn(1); n = mn(2);


%% show images, separate channels
fs = 18; margin = 0.01; cfs = 0.075;
close all
for ci = 1:nc
    means = cellfun(@(x) mean(x,'all'),imgs(:,:,ci));
    [~,I] = max(means,[],'all','linear'); [row,col] = ind2sub(size(means),I);
    lim = stitchedlim(imgs{row,col,ci});
    figure('WindowState','maximized')
    for pidx = 1:npos
        cidx = floor((pidx-1)/ppc) + 1;
        ii = layout(pidx,1); jj = layout(pidx,2);
        if ~isnan(ii)
            idx = sub2ind(size(imgs,[2,1]),jj,ii);
            subplot_tight(size(imgs,1),size(imgs,2),idx,margin)
            imshow(imadjust(imgs{ii,jj,ci},lim),...
                    'InitialMagnification','fit','Border','tight')
            cleanSubplot(fs)
            %channel label
            text(50, m*(1-0.5*cfs), channelLabel{ci},...
                'Color','w','FontUnits','normalized','FontSize',cfs,...
                'FontWeight','bold','Interpreter','none')
            %condition label
            text(50, m*0.5*cfs, conditions{cidx},...
                'Color','w','FontUnits','normalized','FontSize',cfs,...
                'FontWeight','bold','Interpreter','none')
        end
    end
    drawnow
    saveas(gcf,fullfile(dataDir,['overview_',channelLabel{ci},'.png']))
end

%% merged colors
fs = 18; margin = 0.01; cfs = 0.075;
order = [4 2 1 0]; %red, green, blue, gray
colors = {'red','green','blue','white'};

clabel = '';
for ii = 1:length(order)
    if order(ii) > 0
        clabel = strcat(clabel,"\color{",colors{ii},"}",channelLabel{order(ii)}," ");
    end
end

lim = zeros(nc,2);
for ci = 1:nc
    means = cellfun(@(x) mean(x,'all'),imgs(:,:,ci));
    [~,I] = max(means,[],'all','linear'); [row,col] = ind2sub(size(means),I);
    lim(ci,:) = stitchedlim(imgs{row,col,ci});
end

figure('WindowState','maximized')
for pidx = 1:npos
    cidx = floor((pidx-1)/ppc) + 1;
    ii = layout(pidx,1); jj = layout(pidx,2);
    if ~isnan(ii)
        idx = sub2ind(size(imgs,[2,1]),jj,ii);
        subplot_tight(size(imgs,1),size(imgs,2),idx,margin)
        RGB = zeros(m,n,3,'uint16');
        for ri = 1:length(order) %make RGB image
            ci = order(ri);
            if ci > 0
                im = imadjust(imgs{ii,jj,ci},lim(ci,:));
                if ri == 4
                    RGB = RGB + 0.25*repmat(im,1,1,3);
                else
                    RGB(:,:,ri) = RGB(:,:,ri) + 0.75*im;
                end
            end
        end

        imshow(RGB,'InitialMagnification','fit','Border','tight')
        cleanSubplot(fs)

        %condition label
        text(50, m*0.5*cfs, conditions{cidx},...
            'Color','w','FontUnits','normalized','FontSize',cfs,...
            'FontWeight','bold','Interpreter','none')

        %channel label
        if sum(imgs{ii,jj,1},'all')>0
            ypos = m*(1-0.5*cfs); xpos = 0.01*n;
            text(xpos, ypos, clabel,...
                'FontUnits','normalized','FontSize',cfs,'FontWeight','bold')
        end
    end
end

name = 'overview';
for ri = 1:length(order)
    if order(ri) > 0
        name = [name,'_',channelLabel{order(ri)}]; %#ok<AGROW>
    end
end
name = [name,'.png'];
saveas(gcf,fullfile(dataDir,name))
