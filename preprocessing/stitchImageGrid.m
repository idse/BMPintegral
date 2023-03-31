function [stitched, upperleft] = stitchImageGrid(upperleft, imgs)
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
    N = size(imgs{find(~emptycells,1,'first')},1);
    
    totalSize = max(cat(1,upperleft{:})) + N - 1;

%     stitched = zeros(totalSize,'uint16');
    stitched = zeros(totalSize,class(imgs{find(~emptycells,1,'first')}));
    for i = 1:size(imgs,1)
        for j = 1:size(imgs,2)
            if ~isempty(imgs{i,j})
                I = upperleft{i,j}(1):upperleft{i,j}(1)+N-1;
                J = upperleft{i,j}(2):upperleft{i,j}(2)+N-1;
                stitched(I,J) = max(stitched(I,J), imgs{i,j});
            else
                disp(['skipping ' num2str([i j])]);
            end
        end
    end
end