function upperleft = stitchedPreviews(dataDir, meta, options)

if ~isfield(options, 'saveStitchFullSize')
    options.saveStitchFullSize = false;
end
if ~isfield(options, 'type')
	options.type = 'MIP';
end
if ~isfield(options, 'margin')
    margin = 0;
else
    margin = options.margin;
end
if ~isfield(options, 'montages')
    options.montages = 1:meta.nWells;
end

% subsampling factor
if ~isfield(options, 'ss')
    if meta.xSize <= 1024
        ss = 2;
    else
        ss = 4;
    end
else
    ss = options.ss;
end
fullmargin = margin;
margin = round(margin/ss);

MIPfiles = dir(fullfile(dataDir,'MIP',['*' options.type '_*tif']));
s = strsplit(MIPfiles(1).name,['_' options.type]);
barefname = s{1};

gridSize = meta.montageGridSize;
pixelOverlap = round(meta.xSize*meta.montageOverlap/100);
posPerCondition = meta.posPerCondition;

upperleft = {};

for montagenr = options.montages
    
    conditionPositions = posPerCondition*(montagenr-1)+1:posPerCondition*montagenr;
    if isempty(gridSize) 
        gridSize = [round(posPerCondition/2) 2];
    end
    
    %nChannels = numel(MIPfiles)/(meta.nWells*posPerCondition);
    nChannels = meta.nChannels;
    
    for ci = 1:nChannels
        
        disp(['-------------processing channel ' num2str(ci) '--------']);
        
        imgs = {};
        tmax = meta.nTime;
        
        for ti = 1:tmax
            
            disp(['processing time ' num2str(ti)]);

            for pi = conditionPositions

                disp(['reading MIP ' num2str(pi)]);
                % gridSize 1 and 2 may be swapped, I have no way of knowing right now
                [i,j] = ind2sub(gridSize, pi - conditionPositions(1) + 1);

                fname = fullfile(dataDir,'MIP',[barefname sprintf(['_' options.type '_p%.4d_w%.4d.tif'],pi-1,ci-1)]);
                if exist(fname,'file') % to deal with missing files, just add black image
                    imgs{j,i} = imread(fname,ti);
                else
                    warning(['channel ' num2str(ci) ' has MIP missing for position ' num2str(pi) ' (did you make MIP for all channels?)']);
                    imgs{j,i} = zeros([1024 1024],'uint16'); 
                end
            end

            % stitch together
            if ci == 1
                if ~isempty(pixelOverlap) %% && ti == 1
                    % get register positions of upper left corner
                    upperleft{montagenr}{ti} = registerImageGridOld(imgs, pixelOverlap);
                else %if ti == 1 && isempty(pixelOverlap)
                    for pi = conditionPositions
                        [i,j] = ind2sub(gridSize,pi - conditionPositions(1) + 1);
                        upperleft{montagenr}{ti}{j,i} = [1+(j-1)*(meta.ySize + 50), 1+(i-1)*(meta.xSize + 50)];
                    end
                end
            end
            [stitched, upperleft{montagenr}{ti}] = stitchImageGrid(upperleft{montagenr}{ti}, imgs);
            % stitchImageGrid shift upperleft so all images are completely
            % within the stitched image 
            % CAUTION : upperleft = [y x]

            % make clean preview (not for quantitative analysis)
            small = imfilter(stitched,ones(ss)/ss^2);
            small = small(1:ss:end,1:ss:end);
            if ti == 1
                preview = zeros([size(small) tmax],'uint16');
            end
            preview(1:size(small,1), 1:size(small,2), ti) = small;
            
            % save fullsize
            if options.saveStitchFullSize
                
                if ti==1 
                    fullsize = size(stitched);
                end
                
                % set up tag struct for tiff/btf
                tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
                tagstruct.Compression = Tiff.Compression.None;
                tagstruct.BitsPerSample = 16;
                tagstruct.SamplesPerPixel = 1;
                tagstruct.ImageLength = fullsize(1);
                tagstruct.ImageWidth = fullsize(2);
                tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;

                fname = sprintf('stitched_MIP_p%.4d_w%.4d.tif',montagenr-1,ci-1);
                if tmax == 1
                    mode = 'w'; 
                else
                    mode = 'w8'; % big tiff
                end
                % append mode
                if ti > 1
                    mode = 'a';
                end
                    
                t = Tiff(fullfile(dataDir,'MIP',fname), mode);
                
                tmp = zeros(fullsize,'uint16');
                idx = min(fullsize, size(stitched));
                tmp(1:idx(1),1:idx(2)) = stitched(1:idx(1),1:idx(2));
                tmp = tmp(1+fullmargin:end-fullmargin, 1+fullmargin:end-fullmargin);
                if ti > 1
                    t.writeDirectory();
                    t.setTag(tagstruct);
                    t.write(tmp);
                    %imwrite(tmp, fullfile(dataDir,'MIP', fname),'WriteMode','Append');
                else
                    t.setTag(tagstruct);
                    t.write(tmp);
                    %imwrite(tmp, fullfile(dataDir,'MIP', fname));
                end
                t.close();
            end
        end
        
        % set lookup table
        MIPinTime = max(preview,[],3);              
        Ilim = stretchlim(MIPinTime);
        % scale lim back to 16 bit range
        Imin = double(min(MIPinTime(MIPinTime>0))); 
        Imax = round(Ilim(2)*(2^16-1));
        preview = mat2gray(preview, [Imin Imax]);
        preview = uint16((2^16-1)*preview);
        
        % crop
        preview = preview(1+margin:end-margin, 1+margin:end-margin,:); 
            
        bareprevfname = sprintf('stichedPreview_p%.4d_w%.4d',montagenr-1,ci-1);

        % save initial and final frame
        imwrite(mat2gray(preview(:,:,1)),fullfile(dataDir, [bareprevfname '_initial.jpg']));
        if tmax > 1
            
            imwrite(mat2gray(preview(:,:,tmax)),fullfile(dataDir, [bareprevfname '_final.jpg']));
        
            % save as compressed video
            v = VideoWriter(fullfile(dataDir,[bareprevfname '.avi']));%,'Uncompressed AVI');
            v.FrameRate = 5;
            open(v)
            for ti = 1:tmax
                writeVideo(v,mat2gray(preview(:,:,ti)))
            end
            close(v);
        end

%         fname = fullfile(dataDir, [bareprevfname '.tif']);
%         imwrite(preview(:,:,1), fname);
%         for ti = 2:tmax
%             imwrite(preview(:,:,ti), fname,'WriteMode','Append');
%         end
    end
end
end