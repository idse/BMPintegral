function [upperleft,links] = registerImageGrid_v2(imgs, pixelOverlap)
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
    
    
    
    % ---------------------
    % Idse Heemskerk, 2016
    % v2, Seth Teague 2020: solve optimization problem to use redundant
    % shifts to improve robustness
    % ---------------------
    
    % I will assume NxN square images
    % however, not all positions on grid have to contain an image
    % so here I find the first actual image on the grid
    emptycells = cellfun(@isempty,imgs,'UniformOutput',false);
    emptycells = cat(1,emptycells{:});
    N = size(imgs{find(~emptycells,1,'first')},1);
    Np = N - pixelOverlap;
    max_var = 10;
    
    % make upperleft for square grid of images (no stitching)
    grid = [size(imgs,1),size(imgs,2)];
    gridUpperleft = cell(size(imgs));
    for i = 1:size(imgs, 1)
        for j = 1:size(imgs, 2)
            gridUpperleft{i,j} = [1+(i-1)*(N + 50), 1+(j-1)*(N + 50)];
        end
    end
    
    %Make the A matrix for least squares (b ~= Ay)
    n = grid(1);
    D = zeros(n-1,n);
    for j = 1:n-1
        D(j,j) = -1;
        D(j,j+1) = 1;
    end
    A = [kron(eye(n),D);kron(D,eye(n))];
    Ax = A; Ay = A;
    
    %Initialize the b matrices for least squares
    blen = 2*grid(1)*grid(2) - grid(1) - grid(2);
    %for square grid, blen = 2*grid(1)^2 - 2*grid(1)
    b1_x = zeros(blen/2,1);
    b2_x = zeros(blen/2,1);
    b1_y = b1_x;
    b2_y = b2_x;
    %initialize c vectors to hold confidences for displacements
    c1 = b1_x; c2 = b2_x;
    
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
                idx = (grid(2)-1)*(j-1) + i;
                belowshift{i+1,j} = [shifti,shiftj];
                b2_x(idx) = -shiftj;
                b2_y(idx) = -shifti;
                
                % crel is a measure for confidence in the displacement
                if abs(shifti) > 10 || abs(shiftj) > 10
                    crel = 0;
                else
                    crel = cmax/mean(c,'all');
                end
                c2(idx) = crel;
                links(2*i,2*j-1) = crel;
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
                idx = grid(2)*(j-1) + i;
                rightshift{i,j+1} = [shifti,shiftj];
                b1_x(idx) = -shiftj;
                b1_y(idx) = -shifti;
                
                % crel is a measure for confidence in the displacement
                if abs(shifti) > 10 || abs(shiftj) > 10
                    crel = 0;
                else
                    crel = cmax/mean(c,'all');
                end
                c1(idx) = crel;
                links(2*i-1,2*j) = crel;
            end
        end
    end
    
%     cx0 = find(abs([b1_x;b2_x]) > max_var);
%     cy0 = find(abs([b1_y;b2_y]) > max_var);
%     
%     bx = [b1_x + pixelOverlap; b2_x];
%     by = [b1_y; b2_y + pixelOverlap];
%     
%     Ax(cx0,:) = 0;
%     Ay(cy0,:) = 0;
%     bx(cx0) = 0;
%     by(cy0) = 0;

    %vertically concatenate into a single vector
    c = [c1;c2];
    bx = [b1_x + pixelOverlap; b2_x];
    by = [b1_y; b2_y + pixelOverlap];
    %apply confidences/weights
    Ax = Ax.*c;
    Ay = Ay.*c;
    bx = bx.*c;
    by = by.*c;
    
    %calculate the least-squares solution with the Moore-Penrose pseudoinverse
    xposes = pinv(Ax)*bx;
    yposes = pinv(Ay)*by;
    %define the position of the top left image to be (0,0)
    xposes = xposes - xposes(1);
    yposes = yposes - yposes(1);
    %reshape positions into a grid
    xpositions = -transpose(reshape(xposes,grid));
    ypositions = -transpose(reshape(yposes,grid));
    %fix the positions to give the positions of the upper left corners of the
    %images
    for j = 1:size(xpositions,2)-1
        xpositions(:,j+1) = xpositions(:,j+1) + j*1024;
    end
    for j = 1:size(ypositions,1)-1
        ypositions(j+1,:) = ypositions(j+1,:) + j*1024;
    end
    
    %convert positions to integer indices
    xpositions = round(xpositions);
    ypositions = round(ypositions);
    
    %Make upperleft with y (vertical) and x (horizontal) positions
    upperleft = cell(grid);
    for i = 1:grid(1)
        for j = 1:grid(2)
            upperleft{i,j} = [ypositions(i,j), xpositions(i,j)];
        end
    end
