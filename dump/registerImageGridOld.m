function [upperleft, links] = registerImageGridOld(imgs, pixelOverlap)
    % register a grid of overlapping images
    % 
    % [upperleft,links] = registerImageGrid(imgs, pixelOverlap)
    %
    %
    % upperleft:    cell array of positions of upperleft corners with
    %               upperleft corner of upperleft image being (1,1)
    % links:        a matrix containing correlation links between images
    %               empty positions will correspond to cells themselves
    %
    % imgs:         cell array of images 
    %               with rows and cols corresponding to grid
    % pixelOverlap: width of overlapping strip in pixels
    %               if left empty, upperleft for a square grid of images
    %               will be returned
    

    % note: could be made fancier by combining redundant shift information
    
    % ---------------------
    % Idse Heemskerk, 2016
    % ---------------------
    
    % I will assume NxN square images
    % however, not all positions on grid have to contain an image
    % so here I find the first actual image on the grid
    emptycells = cellfun(@isempty,imgs,'UniformOutput',false);
    emptycells = cat(1,emptycells{:});
    N = size(imgs{find(~emptycells,1,'first')},1);
    Np = N - pixelOverlap;
    
    % make upperleft for square grid of images (no stitching)
    gridUpperleft = {};
    for i = 1:size(imgs, 1)
        for j = 1:size(imgs, 2)
            gridUpperleft{i,j} = [1+(i-1)*(N + 50), 1+(j-1)*(N + 50)];
        end
    end
    
    % a matrix containing correlation links between images
    % empty positions will correspond to cells themselves, so
    % 0   1.1  0
    % 1.2  0   1
    % 0   1.3  0 
    % would be a 2x2 grid with relative correlation 1.1 for leftshift
    % between image (1,1) and (1,2), 1.2 for belowshift between (1,1) and
    % (2,1), etc
    links = zeros(2*size(imgs)-1);
    
    % return grid upperleft if there is no overlap to stitch with
    if isempty(pixelOverlap) 
        upperleft = gridUpperleft;
        return;
    end
    
    % otherwise, proceed with calculating displacements in overlap regions
    belowshift = {};
    for i = 1:size(imgs, 1)-1     
        for j = 1:size(imgs, 2)
            
            if ~isempty(imgs{i,j}) && ~isempty(imgs{i+1,j})
                
                img = imgs{i,j}(end-pixelOverlap+1:end,:);
                below = imgs{i+1,j}(1:pixelOverlap,:);
                [shifti,shiftj, cmax, c] = xcorr2fft(below, img);
                belowshift{i+1,j} = [shifti,shiftj];
                
                crel = cmax/mean(c(:));
                links(2*i,2*j-1) = crel;
                disp([num2str(i+1) ' ' num2str(j) ': ' num2str([shifti shiftj]) '   ' sprintf('%0.2e, %.2f',cmax,crel)]);
            end
        end
    end

    disp('right');
    rightshift = {};
    for j = 1:size(imgs, 2)-1
        for i = 1:size(imgs, 1)
            
            if ~isempty(imgs{i,j}) && ~isempty(imgs{i,j+1})

                img = imgs{i,j}(:,end-pixelOverlap+1:end);
                right = imgs{i,j+1}(:,1:pixelOverlap);
                [shifti,shiftj, cmax, c] = xcorr2fft(right, img);
                rightshift{i,j+1} = [shifti,shiftj];
                
                crel = cmax/mean(c(:));
                links(2*i-1,2*j) = crel;
                disp([num2str(i) ' ' num2str(j+1) ': ' num2str([shifti shiftj]) '   ' sprintf('%0.2e, %.2f',cmax,crel)]);
            end
        end
    end

    upperleft = {};
    upperleft{1,1} = [1 1];

    % register each row in the first column
    j = 1;
    for i = 2:size(imgs, 1)
        shift = belowshift{i,j};
        upperleft{i,j} = upperleft{i-1,j} + [Np + shift(1), shift(2)];
    end

    % for each row, register all the columns relative to the first one
    for i = 1:size(imgs, 1)
        for j = 2:size(imgs, 2)
            shift = rightshift{i,j};
            upperleft{i,j} = upperleft{i,j-1} + [shift(1), Np + shift(2)];
        end
    end
end