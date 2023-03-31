function [stitched, upperleft] = stitchImageGridIntensityWeighted(upperleft, imgs)
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
m = size(imgs,1);
n = size(imgs,2);

%summed matrix
S = zeros(totalSize);
%weighted sum matrix
WS = zeros(totalSize);

for ii = 1:m
    for jj = 1:n
        if ~isempty(imgs{ii,jj})
            I = upperleft{ii,jj}(1):upperleft{ii,jj}(1)+N-1;
            J = upperleft{ii,jj}(2):upperleft{ii,jj}(2)+N-1;
            img = double(imgs{ii,jj});
            filtered = medfilt2(img);
            
            WS(I,J) = WS(I,J) + img.*filtered;
            S(I,J) = S(I,J) + filtered;            
        else
            disp(['skipping ' num2str([ii jj])]);
        end
    end
end
% error('check these')

stitched = cast(WS./S,imclass);



end