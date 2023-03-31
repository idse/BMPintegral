function [stitched, upperleft] = stitchImageGridDistWeighted_v2(upperleft, imgs)
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
% Updated to take a weighted average over the overlap region instead of
% taking the maximum - Seth Teague, 2021
% ---------------------

% note: could be made fancier by averaging the overlap

UL = cat(1,upperleft{:});
dUL = 1 - min(UL);
for i = 1:numel(upperleft)
    if ~isempty(upperleft{i})
        upperleft{i} = upperleft{i} + dUL;
    end
end

% again assuming square images
% however, not all positions on grid have to contain an image
% so here I find the first actual image on the grid
emptycells = cellfun(@isempty,imgs,'UniformOutput',false);
emptycells = cat(1,emptycells{:});
imsize = size(imgs{find(~emptycells,1,'first')},[1,2]);
M = imsize(1); N = imsize(2);

totalSize = max(cat(1,upperleft{:})) + [M - 1, N - 1];

m = size(imgs,1);
n = size(imgs,2);

imclass = class(imgs{find(~emptycells,1,'first')});
stitched = zeros(totalSize);

for ii = 1:m
    for jj = 1:n
        upper = upperleft{ii,jj}(1);
        left = upperleft{ii,jj}(2);
        if ii == 1
            boverlap = M - (upperleft{ii+1,jj}(1) - upper); %bottom overlap
            leftweight = diag([ones(1,M-boverlap),linspace(1,0,boverlap)]);
        elseif ii == m
            toverlap = M - (upper - upperleft{ii-1,jj}(1)); %top overlap
            leftweight = diag([linspace(0,1,toverlap),ones(1,M-toverlap)]);
        else
            boverlap = M - (upperleft{ii+1,jj}(1) - upper);
            toverlap = M - (upper - upperleft{ii-1,jj}(1));
            leftweight = diag([linspace(0,1,toverlap),...
                ones(1,M-toverlap-boverlap),linspace(1,0,boverlap)]);   
        end

        if jj == 1
            roverlap = N - (upperleft{ii,jj+1}(2) - left); %right overlap
            rightweight = diag([ones(1,N-roverlap),linspace(1,0,roverlap)]);
        elseif jj == n
            loverlap = N - (left - upperleft{ii,jj-1}(2)); %left overlap
            rightweight = diag([linspace(0,1,loverlap),ones(1,N-loverlap)]);
        else
            roverlap = N - (upperleft{ii,jj+1}(2) - left); %right overlap
            loverlap = N - (left - upperleft{ii,jj-1}(2)); %left overlap
            rightweight = diag([linspace(0,1,loverlap),...
                ones(1,N-loverlap-roverlap),linspace(1,0,roverlap)]);
        end
        
        I = upperleft{ii,jj}(1):upperleft{ii,jj}(1)+M-1;
        J = upperleft{ii,jj}(2):upperleft{ii,jj}(2)+N-1;
        
        stitched(I,J) = stitched(I,J) + leftweight*double(imgs{ii,jj})*rightweight;
        
    end
end

stitched = cast(stitched,imclass);


end