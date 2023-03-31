function Ilim = makeXsection(im, fname, meta, index, RGBset, nucChannel, tol, zcorrection, Ilim)
    % im: image stack
    % meta: metadata
    % index : y index to make cross section at
    % RGBset : channels to use for RGB image

    lw = 10;

    if ~exist('zcorrection','var')
        zcorrection = false;
    end
    
    % ------- z-correction
    if zcorrection
        % sub bg
        imbgsub = im;
        tic
        for zi = 1:size(im,4)
            for ci = 1:4
                slice = imgaussfilt(imbgsub(:,:,ci,zi), 2);
                bg = imopen(slice, strel('disk',50));
                imbgsub(:,:,ci,zi) = slice - bg;
            end
        end
        toc

        % correct intensity
        imcorr = imbgsub;
        f = [1.1 1.1 1.1, 1.1];
        tic
        for zi = 1:size(imcorr,4)
            for ci = 1:4
                imcorr(:,:,ci,zi) = imcorr(:,:,ci,zi)*(f(ci)^(zi-1));
            end
        end
        im = imcorr;
        toc
    end
    % ------- 
    
    if ~exist('Ilim','var')
        Ilim = {};
        fillIlim = true;
    else
        fillIlim = false;
    end
    MIP = max(im,[],4);

    MIPadj = MIP; 
    for ci = 1:meta.nChannels
        if fillIlim
            Ilim{ci} = stitchedlim(MIP(:,:,ci), struct('Tol',tol(ci,:)));
        end
        MIPadj(:,:,ci) = imadjust(MIP(:,:,ci),Ilim{ci}); 
    end

    t = round(lw/2);
    MIPadj(index-t:index+t,:,:) = max(MIPadj(:));

    newnslices = round(meta.zres / meta.xres)*size(im,4);
    scaled = zeros([newnslices size(im,2)],'uint16');
    for ci = 1:meta.nChannels
        scaled(:,:,ci) = flipud(imresize(squeeze(max(im(index,:,ci,:),[],1)), [size(im,2) newnslices])');
        scaled(:,:,ci) = imadjust(scaled(:,:,ci),Ilim{ci});
    end

    %----- display

    figure,
    xyim = MIPadj(:,:,RGBset);
    if ~isempty(nucChannel)
        xyimnuc = MIPadj(:,:,nucChannel);
        for i = 1:3
           xyim(:,:,i) =  xyim(:,:,i) + xyimnuc;
        end
    end
    imshow(xyim,[])
    imwrite(xyim, [fname(1:end-4) '_' meta.channelLabel{RGBset} '_' num2str(index([1 end])) '_MIP.png']);

    figure,
    for ci = 1:meta.nChannels

        subplot_tight(meta.nChannels + 1, 1, ci)
        imshow(scaled(:,:,ci),[])
        title(meta.channelLabel{ci})
    end

    subplot_tight(meta.nChannels + 1, 1, meta.nChannels + 1)
    xzim = scaled(:,:,RGBset);
    if ~isempty(nucChannel)
        xzimnuc = scaled(:,:,nucChannel);
        for i = 1:3
           xzim(:,:,i) =  xzim(:,:,i) + xzimnuc;
        end
    end
    imshow(xzim,[])

    str = [];
    for i = 1:3
        str = [str meta.channelLabel{RGBset(i)} '; '];
    end
    title(str);
    saveas(gcf, [fname(1:end-4) '_xsect_' num2str(index([1 end])) '.png']);
end