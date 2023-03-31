function Ilim = micropatternPieVis(dataDir, position, options)

    if ~isfield(options,'tol')
        tol = cat(2, [1 1 1 1]'*0.01, [1 1 1 1]'*0.99);
    else
        tol = options.tol;
    end
    if ~isfield(options,'channels')
        channels = 2:4;
    else
        channels = options.channels;
    end
    if ~isfield(options,'pieOrder')
        order = 1:3;
    else
        order = options.pieOrder;
    end
    if ~isfield(options,'scalebar')
        options.scalebar = true;
    end
    if ~isfield(options,'segOverlay')
        options.segOverlay = false;
    end
    if isfield(options,'margin')
        cropmargin = options.margin;
    else
        cropmargin = 150; 
    end
    
    margin = 300;
    PI = atan(1)*4;

    s = strsplit(position.filename,{'_FusionStitcher','_Stitched','.ims','.tif'});
    prefix = [s{:}];
    subDir = [prefix '_zslices'];

    imgs = {};
    Ilim = {};
    for cii = 1:4

        img = imread(fullfile(dataDir, subDir, sprintf([prefix '_MIP_w%.4d.tif'], cii-1)));
        
        if ~isfield(options,'Ilim')
            Ilim{cii} = stretchlim(img,tol(cii,:));
        else
            Ilim{cii} = options.Ilim{cii};
        end
        imgs{cii} = imadjust(img,Ilim{cii});

        %add margin
        X = zeros(size(img)+2*margin);
        X(margin+1:margin+size(img,1),margin+1:margin+size(img,2)) = imgs{cii};
        imgs{cii} = X;
    end

    xres = position.radiusMicron / position.radiusPixel;
    center = position.center + margin;
    Rmax = uint16(position.radiusPixel + 25/xres); %25 micron margin
    [X,Y] = meshgrid(1:size(imgs{1},2),1:size(imgs{1},1));
    R = sqrt((X - center(1)).^2 + (Y - center(2)).^2);
    F = atan2((X - center(1)),-(Y - center(2)));
    disk = R > Rmax;

    imgs_tmp = imgs;
    for cii = 1:4
        imgs_tmp{cii} = mat2gray(imgs_tmp{cii});
        imgs_tmp{cii}(disk) = 1; % change 1 to 0 to get black background
    end
    
    % centered crop indices
    CM = uint16(center);
    Rcrop = Rmax + round(cropmargin/2);
    yrange = CM(2)-Rcrop:CM(2)+Rcrop;
    xrange = CM(1)-Rcrop:CM(1)+Rcrop;
    
    % 50um scalebar on DAPI image
    Npixels = uint16(50/xres);
    marg = 50;
    w = 8;
    DAPIim = mat2gray(imgs_tmp{1}(yrange,xrange));
    if options.scalebar
        DAPIim(end-marg:end-marg+w, marg:marg+Npixels) = 0;
    end
    imwrite(DAPIim, fullfile(dataDir, sprintf([prefix '_DAPI.png'])));
        
    %-----------------------
    if options.segOverlay
        fname = sprintf([position.filename(1:end-4) '_SegOverlay.png']);
        segOverlay_in = imread(fullfile(dataDir, fname));
        X = zeros([size(img)+2*margin 3],'uint8');
        X(margin+1:margin+size(img,1),margin+1:margin+size(img,2),:) = segOverlay_in;
        segOverlay_out = X;
        for i = 1:3
            im = segOverlay_out(:,:,i);
            im(disk)=255;
            segOverlay_out(:,:,i) = im;
        end
        imwrite(segOverlay_out(yrange,xrange,:), fullfile(dataDir, sprintf([prefix '_SegOverlay_centered.png'])));
    end
    %-----------------------
    
    imgs_tmp = imgs_tmp(channels);
    overlay = cat(3,imgs_tmp{:});
    overlay = overlay(yrange,xrange,:);
    if options.scalebar
        overlay(end-marg:end-marg+w, marg:marg+Npixels,:) = 0;
    end
    imwrite(overlay, fullfile(dataDir, sprintf([prefix '_combinedRGB.png'])));
    
    for ci = 1:3
        imwrite(imgs_tmp{ci}(yrange,xrange), fullfile(dataDir, sprintf([prefix '_w%.4d.png'],ci)));
    end
    
    %-----------------------
    mask = {}; 
    N = 3;
    for i = 1:N
        mask{i} = F < -PI + (i-1)*2*PI/N | F > -PI + i*2*PI/N;
        mask{i} = imdilate(mask{i},strel('disk',10));
    end
    lines = mask{1};
    for i = 2:N
        lines = lines & mask{i};
    end

    imgs_tmp = imgs(channels);
    imgs_tmp = imgs_tmp(order);
    for cii = 1:numel(channels)
        imgs_tmp{cii} = mat2gray(imgs_tmp{cii});
        imgs_tmp{cii}(mask{cii}) = 0;
    end
    combined = sum(cat(3,imgs_tmp{:}),3);
    combined(disk) = 1;
    combined(lines) = 1;
    combined = combined(yrange,xrange,:);
    %imshow(combined,[])
    if options.scalebar
        combined(end-marg:end-marg+w, marg:marg+Npixels) = 0;
    end
    imwrite(combined, fullfile(dataDir, sprintf([prefix '_combined.png'])));

    % % color version
    % 
    % imgs_tmp = imgs;
    % for cii = 1:numel(channels)
    %     imgs_tmp{cii} = mat2gray(imgs_tmp{cii});
    %     imgs_tmp{cii}(mask{cii}) = 0;
    %     imgs_tmp{cii}(disk) = 1;
    %     imgs_tmp{cii}(lines) = 1;
    % end
    % combined = cat(3,imgs_tmp{:});
    % imshow(combined,[])
    
end