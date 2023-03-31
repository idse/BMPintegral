function makeSegmentationCrossSection(img,LM,meta,yinds,xinds)

%img = nuclear channel image stack
%L = 3D label matrix
%meta = experiment metadata
%xinds = x indices
%yinds = y indices

if numel(xinds) > 1 && numel(yinds) > 1
    error('either x or y range must be only a single index')
end

if numel(xinds) == 1
    res = meta.yres;
else
    res = meta.xres;
end

im = imadjust(squeeze(img(xinds,yinds,:))');
L = squeeze(LM(xinds,yinds,:))';
im2 = visualize_nuclei_v2(L,im);

ratio = (size(im,1)*meta.zres)/(size(im,2)*res);
newsize = round(2048*[ratio,1]);

im = flipud(imresize(im,newsize));
im2 = flipud(imresize(im2,newsize));

figure
subplot_tight(2,1,1)
imshow(im)
title('Image')
%cleanSubplot

subplot_tight(2,1,2)
imshow(im2)
title('With Labels')
%cleanSubplot

end