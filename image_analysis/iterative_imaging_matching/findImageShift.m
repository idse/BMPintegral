function shiftyx = findImageShift(im1, im2, alignmentMethod)
%pi is the index of the position for which the images are being aligned
%this is mostly called from within mapPoints instead of as a
%standalone function
%XY is a ncellsx2 array of adjusted cell locations in fixed image

if ~exist('alignmentMethod','var')
    alignmentMethod = 'automatic';
end

im1 = max(im1,[],3); %take MIP in case of z-slices
im2 = max(im2,[],3); %take MIP in case of z-slices
%make images the same size
mn = min(size(im1), size(im2));
m = mn(1); n= mn(2);
I1 = im1(1:m,1:n);
I2 = im2(1:m,1:n);

if strcmp(alignmentMethod,'automatic')
    %automatically determine shift between fixed and live images
    [shiftx,shifty,~,~] = xcorr2fft(I1, I2);
    shiftyx = [shifty shiftx];
elseif strcmp(alignmentMethod,'manual')
    I1 = imadjust(I1,stitchedlim(I1));
    I2 = imadjust(I2,stitchedlim(I2));
    shiftx = 0;
    shifty = 0;
    close all
    f = figure('WindowState','maximized');
    p1 = imshow(cat(3,I1,I2,I1),'InitialMagnification','fit');
    cleanSubplot
    breakvar = false;
    while breakvar == false
        xinrange = max(1,1+shiftx):min(n,n+shiftx);
        xoutrange = max(1,1-shiftx):min(n,n-shiftx);
        yinrange = max(1,1+shifty):min(m,m+shifty);
        youtrange = max(1,1-shifty):min(m,m-shifty);

        aligned = zeros(mn,'uint16');
%         aligned(xoutrange, youtrange) = I1(xinrange, yinrange);
        aligned(youtrange, xoutrange) = I1(yinrange, xinrange);
        %update image for new frame
        set(p1,'CData',cat(3,aligned,I2,aligned));

        title(sprintf('shiftx = %d, shifty = %d', shiftx, shifty))

        waitforbuttonpress;
        key = f.CurrentCharacter;

        if strcmp(key,'d')
            shiftx = shiftx - 1;
        elseif strcmp(key,'a')
            shiftx = shiftx + 1;
        elseif strcmp(key,'s')
            shifty = shifty - 1;
        elseif strcmp(key,'w')
            shifty = shifty + 1;
        elseif strcmp(key,'+')
            xl = get(gca,'xlim');
            yl = get(gca,'ylim');
            set(gca,'xlim', mean(xl) + 0.25*(xl(2) - xl(1))*[-1 1])
            set(gca,'ylim', mean(yl) + 0.25*(yl(2) - yl(1))*[-1 1])
        elseif strcmp(key,'-')
            xl = get(gca,'xlim');
            yl = get(gca,'ylim');
            set(gca,'xlim', mean(xl) + (xl(2) - xl(1))*[-1 1])
            set(gca,'ylim', mean(yl) + (yl(2) - yl(1))*[-1 1])
        elseif strcmp(key,'8')
            shifty = shifty + 10;
        elseif strcmp(key,'2')
            shifty = shifty - 10;
        elseif strcmp(key,'6')
            shiftx = shiftx - 10;
        elseif strcmp(key,'4')
            shiftx = shiftx + 10;
        elseif strcmp(key,'e')
            breakvar = true;
        end
        drawnow limitrate
    end
    close(f)
    shiftyx = [-shiftx, -shifty];
end


end