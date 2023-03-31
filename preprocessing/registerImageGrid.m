function [upperleft,links] = registerImageGrid(imgs, pixelOverlap)
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
                
                % crel is a measure for confidence in the displacement
                crel = cmax/mean(c(:));
                links(2*i,2*j-1) = crel;
                %disp([num2str(i+1) ' ' num2str(j) ': ' num2str([shifti shiftj]) '   ' sprintf('%0.2e, %.2f',cmax,crel)]);
            end
        end
    end

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
                %disp([num2str(i) ' ' num2str(j+1) ': ' num2str([shifti shiftj]) '   ' sprintf('%0.2e, %.2f',cmax,crel)]);
            end
        end
    end

	% registerImageGridOld works well for square grid, 
    % but not if you leave images out
    % code below is harder to follow but works when corner images are
    % missing (are empty for high-res micropattern data)
    upperleft = {};
    
    % register rows in each column
    for j = 1:size(imgs, 2)
        
        shift = [0 0];
        for i = 1:size(imgs, 1)

            % initialize 
            upperleft{i,j} = [(i-1)*Np + 1, (j-1)*Np + 1];

            % then update
            if ~isempty(belowshift{i,j})
                shift = shift + belowshift{i,j};
                upperleft{i,j} = upperleft{i,j} + shift;
            end
        end
    end
    
    % register the columns
    shift = [0 0];
    for j = 2:size(imgs, 2)

        % find the first image in the column that is to the right of
        % another image, and use that for the shift
        i = 1;
        while isempty(rightshift{i,j}) i = i+1; end
        shift = shift + rightshift{i,j};

        for i = 1:size(imgs,1)
            upperleft{i,j} = upperleft{i,j} + shift;
        end
    end
end