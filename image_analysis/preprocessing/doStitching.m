clear; clc; close all

%% setup
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
% dataDir = 'Y:\Seth\220315_Rues2Smad1RFP_BMP_IWP2_live\nikon';
dataDir = scriptPath;
bare = '_p%.4d_w%.4d_t%.4d.tif';

mipdir = fullfile(dataDir, 'MIP');
if ~exist(mipdir,'dir'), mkdir(mipdir); end

gridsize = [2 3];
ppc = gridsize(1)*gridsize(2); %positions per condition (stitched grid)
montageOverlap = 20;
nucChannel = 0;
montageFusionGrid = false;

%possible extensions: .nd2, .ims, .lif, .tif
ext = '.nd2';

listing = dir(fullfile(dataDir,['*',ext]));
fname = fullfile(dataDir,listing(1).name);
r = bfGetReader(fname);
m = r.getSizeY; n = r.getSizeX; nt = r.getSizeT; nz = r.getSizeZ; nc = r.getSizeC;
pixelOverlapY = round(m*montageOverlap/100); pixelOverlapX = round(n*montageOverlap/100);

stitchedSize = round([gridsize(1) - montageOverlap*(gridsize(1)-1)/100,... 
                        gridsize(2) - montageOverlap*(gridsize(2)-1)/100].*[m,n]);

if strcmp(ext,'.nd2') || strcmp(ext,'.lif')
    %do we also need to account for multiple files with multiple series?
    npos = r.getSeriesCount;
    seriesflag = true;
else
    npos = length(listing);
    seriesflag = false;
end

ngrids = npos/ppc;

%populate and save metadata
% manualMeta = struct(...
%     'nWells',   ngrids,...
%     'posPerCondition',  ppc,...
%     'montageGridSize',  gridsize,...
%     'montageOverlap',   montageOverlap);
% meta = Metadata(dataDir, manualMeta);
% save(fullfile(dataDir,'meta.mat'),'meta');

fprintf('xyzct = %d,%d,%d,%d,%d\n',n,m,nz,nc,nt)
fprintf('number of positions = %d\n',npos)
fprintf('number of grids = %d\n',ngrids)

%determine order of images in a stitched grid
order = NaN(ppc,2);
idx = 1;
if ~montageFusionGrid
    for ii = 1:gridsize(1)
        if mod(ii,2) == 0 %alternate direction for odd and even rows
            roworder = gridsize(2):-1:1;
        else
            roworder = 1:gridsize(2);
        end
        order(idx:idx + gridsize(2) - 1,:) = [ii*ones(gridsize(2),1),roworder(:)];
        idx = idx + gridsize(2);
    end
else
    for ii = 1:gridsize(1)
        roworder = 1:gridsize(2);
        order(idx:idx + gridsize(2) - 1,:) = [ii*ones(gridsize(2),1),roworder(:)];
        idx = idx + gridsize(2);
    end
end


grid = NaN(gridsize);
disp(order)
for idx = 1:ppc
    ii = order(idx,1); jj = order(idx,2);
    grid(ii,jj) = idx;
end
disp(grid)

%% do stitching

tic
for pidx = 1:ngrids
    fprintf('Stitching grid #%d\n',pidx)
    for ti = 1:nt
        fprintf('t = %d\n',ti)
        imgs = cell(gridsize);
        nucmips = cell(gridsize);
        for idx = 1:ppc
            %overall index of the (idx)th image in the (pidx)th grid
            fidx = (pidx - 1)*ppc + idx;
            fprintf('.')
            if seriesflag
                %if all positions are saved as series in a single file
                r.setSeries(fidx-1);
            else
                %if each position is saved in a separate file
                %it would be faster not to have to reinitialize these
                %readers at each time point, but it would take up too much
                %ram to load all positions into a 
                fname = fullfile(dataDir,listing(fidx).name);
                r = bfGetReader(fname);
            end

            ii = order(idx,1); jj = order(idx,2);
            img = zeros(m,n,nz,nc,'uint16');
            for ci = 1:nc
                for zi = 1:nz
                    im = bfGetPlane(r, r.getIndex(zi-1,ci-1,ti-1)+1);
                    img(:,:,zi,ci) = im;
                end
            end
            imgs{ii,jj} = img;
            nucmips{ii,jj} = max(squeeze(img(:,:,:,nucChannel+1)),[],3);
        end
        fprintf('\n')

        if ti == 1 %double check the overlap values in pixels on the first time point
            %check pixelOverlapY
            [~,I] = max((cellfun(@(x) mean(x(end-pixelOverlapY+1:end,:),'all'),nucmips(1:end-1,:))),[],'all','linear');
            [ii,jj] = ind2sub(gridsize - [1 0],I);
            im1 = nucmips{ii,jj}(end-pixelOverlapY+1:end,:);
            im2 = nucmips{ii+1,jj}(1:pixelOverlapY,:);
            [offset, ~, ~, ~] = xcorr2fft(im1, im2);

            if abs(offset/pixelOverlapY) > 0.25
                warning('detected pixel overlap very different from nominal value, sticking with nominal value');
            else
                pixelOverlapY = pixelOverlapY + offset;
            end

            %check pixelOverlapX
            [~,I] = max((cellfun(@(x) mean(x(:,end-pixelOverlapX+1:end),'all'),nucmips(:,1:end-1))),[],'all','linear');
            [ii,jj] = ind2sub(gridsize - [0 1],I);
            im1 = nucmips{ii,jj}(:,end-pixelOverlapX+1:end);
            im2 = nucmips{ii,jj+1}(:,1:pixelOverlapX);
            [~, offset, ~, ~] = xcorr2fft(im1, im2);

            if abs(offset/pixelOverlapX) > 0.25
                warning('detected pixel overlap very different from nominal value, sticking with nominal value');
            else
                pixelOverlapX = pixelOverlapX + offset;
            end
        end

        upperleft = registerImageGrid_v3(nucmips, [pixelOverlapY pixelOverlapX]);
        for ci = 1:nc
            suffix = sprintf(bare,pidx-1,ci-1,ti-1);
            fname = fullfile(dataDir,['stitched',suffix]);
            mipname = fullfile(mipdir,['stitched_MIP',suffix]);
            idxname = fullfile(mipdir,['stitched_MIPidx',suffix]);

            stitched = zeros(stitchedSize(1),stitchedSize(2),nz,'uint16');
            for zi = 1:nz
                if zi == 1
                    mode = 'overwrite';
                else
                    mode = 'append';
                end
                ims = cellfun(@(x) x(:,:,zi,ci),imgs,'UniformOutput',false);
                [stitch, ~] = stitchImageGridDistWeighted(upperleft, ims);

                %trim to projected size
                tmp = zeros(stitchedSize,'uint16');
                yrange = 1:min(stitchedSize(1),size(stitch,1));
                xrange = 1:min(stitchedSize(2),size(stitch,2));
                tmp(yrange, xrange) = stitch(yrange, xrange);
                stitched(:,:,zi) = tmp;

                imwrite(tmp,fname,'WriteMode',mode)
            end
            [MIP,MIPidx] = max(stitched,[],3);
            imwrite(MIP,mipname)
            if ci == nucChannel + 1
                MIPidx = uint8(MIPidx);
                imwrite(MIPidx,idxname)
            end
        end


    end
end
toc


%% write stitched previews
bare = 'stitched_MIP_p%.4d_w%.4d_t%.4d.tif';
writebare = 'preview_p%.4d_w%.4d.mp4';

for pidx = 1:ngrids
    fprintf('Position %d\n',pidx)
    for ci = 1:nc
        fprintf('Channel %d\n',ci)
        writename = fullfile(dataDir,sprintf(writebare,pidx-1,ci-1));
        v = VideoWriter(writename,'MPEG-4'); %#ok<TNMLP>
        v.FrameRate = 5;
        v.Quality = 100;
        open(v)
        for ti = 1:nt
            fprintf('.')
            if mod(ti,45) == 0
                fprintf('\n')
            end
            name = fullfile(mipdir,sprintf(bare,pidx-1,ci-1,ti-1));
            img = imread(name);
            if ti == 1
                lim = seglim(img);
%                 lim = stitchedlim(img);
            end
            img = imread(name);
            writeVideo(v,mat2gray(imadjust(img,lim)))
        end
        fprintf('\n')
        close(v)
    end
end

