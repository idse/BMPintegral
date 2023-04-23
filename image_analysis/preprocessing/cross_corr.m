function [offset, corrs] = cross_corr(im1,im2,dim,k)
%perform zero-normalized cross correlation
%k is the largest offset checked
%dim is the dimension in which to shift the images
%im1 is the first image; it should be either above or to the left of im2

if ~exist('k','var')
    k = min(size(im1,dim),size(im2,dim));
end

%If input images are in color first convert them to grayscale
if size(im1,3) > 1
    im1 = rgb2gray(im1);
end

if size(im2,3) > 1
    im2 = rgb2gray(im2);
end

%Center the distribution of 
% im1 = im1 - mean(mean(im1));
% im1 = im1/max(max(abs(im1)));
% im2 = im2 - mean(mean(im2));
% im2 = im2/max(max(abs(im2)));

corrs = Inf(k,1); %stores the cross-correlation value for each shift from 1 to k
if dim == 1
    for j=1:k
        %extract the overlapping portion of the images after shifting
        test1 = im1(end-j+1:end,:);
        %convert to a vector and subtract the mean to center the data about
        %zero
        test1 = single(test1(:) - mean(mean(test1)));
        %normalize the vector
        test1 = test1/norm(test1);
        
        test2 = im2(1:j,:);
        test2 = single(test2(:) - mean(mean(test2)));
        test2 = test2/norm(test2);
        
        corrs(j) = sum(test1.*test2,'all');%/numel(test1);
    end
elseif dim == 2
    for j=1:k
        test1 = im1(:,end-j+1:end);
        test1 = single(test1(:) - mean(mean(test1)));
        test1 = test1/norm(test1);
        
        test2 = im2(:,1:j);
        test2 = single(test2(:) - mean(mean(test2)));
        test2 = test2/norm(test2);
        corrs(j) = sum(test1.*test2,'all');%/numel(test1);
    end
end

[~, offset] = max(corrs);

end
