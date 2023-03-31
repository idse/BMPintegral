function segoverlay = nuclearSegmentationOverlay(im, seg, fgseg)
    
    if ~exist('fgseg','var') || isempty('fgseg')
        fgseg = zeros([size(seg) 3]);
    else
        fgseg = 255*(fgseg - seg);
    end
    
    L = bwlabel(seg);
    
    N = max(L(:));
    cmap = jet(N);
    p = randperm(N);
    cmap = cmap(p,:);
    
    Lrgb = label2rgb(L,cmap,'k');
    
    if ~isempty(im)
        I = imadjust(mat2gray(im), stretchlim(mat2gray(im),0.005));
        segoverlay = uint8(double(I)*255*0.5 + 0.4*fgseg + 0.7*double(Lrgb));
    else
        segoverlay = uint8(0.8*fgseg + double(Lrgb));
    end
end