function makeCrossSection(img,LM,meta,yinds,xinds)

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

MIP = max(img,[],3);
figure
imshow(imadjust(MIP,stitchedlim(MIP)));
y = [min(xinds),max(xinds)];
x = [min(yinds),max(yinds)];
line(x,y,'Color','r','LineWidth',2)

im = imadjust(squeeze(img(yinds,xinds,:))');
L = squeeze(LM(yinds,xinds,:))';
im2 = visualize_nuclei_v2(L,im);

ratio = (size(im,1)*meta.zres)/(size(im,2)*res);
newsize = round(2048*[ratio,1]);

im = imresize(im,newsize);
im2 = imresize(im2,newsize);

figure
subplot_tight(2,1,1)
imshow(im)
title('Image')
cleanSubplot

subplot_tight(2,1,2)
imshow(im2)
title('With Labels')
cleanSubplot

end