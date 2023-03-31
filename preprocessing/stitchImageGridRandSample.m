function [stitched, upperleft] = stitchImageGridRandSample(upperleft, imgs)
% combine image grid into single image
%
% stitched = stitchImageGrid(upperleft, imgs)
%
%
% stitched:     combined image
%
% upperleft:    cell array of positions of upperleft corners with
%               upperleft corner of upperleft image being (1,1)
%               as produced by registerImageGrid(..)
% imgs:         cell array of images 
%               with rows and cols corresponding to grid
%
% see also registerImageGrid

% ---------------------
% Idse Heemskerk, 2016
% ---------------------

% note: could be made fancier by averaging the overlap

UL = cat(1,upperleft{:});
dUL = 1 - min(UL);
for ii = 1:numel(upperleft)
    if ~isempty(upperleft{ii})
        upperleft{ii} = upperleft{ii} + dUL;
    end
end

% again assuming square images
% however, not all positions on grid have to contain an image
% so here I find the first actual image on the grid
emptycells = cellfun(@isempty,imgs,'UniformOutput',false);
emptycells = cat(1,emptycells{:});
N = size(imgs{find(~emptycells,1,'first')},1);

totalSize = max(cat(1,upperleft{:})) + N - 1;

imclass = class(imgs{find(~emptycells,1,'first')});
stitched = zeros(totalSize,imclass);
ims = cell(1,1,4);
for ii = 1:4
    ims{ii} = zeros(totalSize,imclass);
end

m = size(imgs,1);
n = size(imgs,2);
Z = zeros(m,n);
starts = [1,1;1,2;2,1;2,2];
for idx = 1:4
    for ii = starts(idx,1):2:m
        for jj = starts(idx,2):2:n
            Z(ii,jj) = idx;
        end
    end
end

for ii = 1:m
    for jj = 1:n
        if ~isempty(imgs{ii,jj})
            I = upperleft{ii,jj}(1):upperleft{ii,jj}(1)+N-1;
            J = upperleft{ii,jj}(2):upperleft{ii,jj}(2)+N-1;
            ims{Z(ii,jj)}(I,J) = imgs{ii,jj};
            stitched(I,J) = imgs{ii,jj};
        else
            disp(['skipping ' num2str([ii jj])]);
        end
    end
end

% M = totalSize(1); N = totalSize(2);
ims = cell2mat(ims);
% error('temperror')
Nsamps = sum(ims>0,3);
for ll = 2:4
    inds = find(Nsamps==ll);
    perminds = randi(ll,length(inds),1);
    for idx = 1:length(inds)
        ind = inds(idx);
        [ii,jj] = ind2sub(totalSize,ind);
        vals = ims(ii,jj,:); vals = vals(vals>0);
        stitched(ii,jj) = vals(perminds(idx));
    end
end




end
