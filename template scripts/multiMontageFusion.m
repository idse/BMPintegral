clear all; close all;
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);

% raw data location, modify if not the same as the location of this script
% protocolFile = '/Users/8R6KGQ2/desktop/210215_LL1_ESI_ECad_0_hr.bkp';
protocolFile = 'Y:\Seth\Fusion protocols\210216_BMP_IWP2_fixed.bkp';

% um/pixel
resolution = 0.301;      %0.603; %20x         %1.26; 10x?           0.301; % 40x
micropattern  = false;

fit_z_to_plane = true; %for uneven plate/grid, select multiple points
pointsPerGrid = 3; %any number you want that is >= 3

gridSize = [6 6]; % for 40x
montageOverlap = 20; % percentage

% 3x3 grid will have length (2*(1024-overlap) + 1024-2*overlap)*resolution
% at 40x that is 868 micron so large enough for a 700 micron colony
overlapPixel = round(1024*montageOverlap/100);

%% make grid centered around single position (or average of two)

protocol = readFusionProtocol(protocolFile);

XYZ = protocol.XYZ;
nGridPositions = gridSize(1)*gridSize(2);

gridXYZ = {};
AF = {};
if micropattern
    nColonies = size(XYZ,1)/2;
elseif fit_z_to_plane
    nColonies = size(XYZ,1)/pointsPerGrid;
    if mod(nColonies,1) > 0
        error('wrong number of positions; must be a multiple of pointsPerGrid')
    end
else
    nColonies = size(XYZ,1);
end

if fit_z_to_plane && pointsPerGrid < 3
    %least squares fit to equation ax + by + c = z to picked positions
    coeffs = pinv([XYZ(:,1),XYZ(:,2),ones(size(XYZ,1),1)])*XYZ(:,3);
    a = coeffs(1); b = coeffs(2); c = coeffs(3);
end

for pi = 1:nColonies

    if micropattern
        x = mean(XYZ(2*pi-1:2*pi,:),1);
        af = round(mean(protocol.AFoffset(2*pi-1:2*pi)));
    elseif fit_z_to_plane
        xyz = XYZ(pointsPerGrid*(pi-1)+1:pointsPerGrid*pi,:);
        x = (min(xyz,[],1) + max(xyz,[],1))/2;
        af = round(mean(protocol.AFoffset(pointsPerGrid*(pi-1)+1:pointsPerGrid*pi)));
    else
        x = XYZ(pi,:);
        af = protocol.AFoffset(pi);
    end
    
    % Fusion specifies xy in mm, not um
    spacing = (1024 - overlapPixel)*resolution/1000; 

    gridXYZ{pi} = zeros([nGridPositions 3]);
    AF{pi} = zeros([nGridPositions 1]);

    nmax = (gridSize(1)-1)/2;
    mmax = (gridSize(2)-1)/2;
    
    i = 1;
    for n = -nmax:nmax
        for m = -mmax:mmax
            gridXYZ{pi}(i,:) = [x(1) + m*spacing, x(2) + n*spacing, x(3)];
            AF{pi}(i) = af;
            i = i+1;
        end
    end
    
    if fit_z_to_plane
        if pointsPerGrid > 2
            %fit the plane only to the points in this well if using 3 or
            %more positions per grid
            coeffs = pinv([xyz(:,1),xyz(:,2),ones(pointsPerGrid,1)])*xyz(:,3);
            a = coeffs(1); b = coeffs(2); c = coeffs(3);
        end
        %x and y of grid positions
        xy = gridXYZ{pi}(:,1:2);
        %use plane to determine z for each position in the grid
        z = a*xy(:,1) + b*xy(:,2) + c;
        gridXYZ{pi}(:,3) = z;
    end
end

gridXYZcombined = cat(1,gridXYZ{:});
AFCombined = cat(1,AF{:});

% visualize
w = 1024*resolution/1000;
h = w;
for i = 1:size(gridXYZcombined,1)
    rectangle('Position',[gridXYZcombined(i,1)-w/2,gridXYZcombined(i,2)-h/2,w,h])
end
hold on
scatter(gridXYZ{1}(:,1), gridXYZ{1}(:,2));
scatter(XYZ(:,1),XYZ(:,2),'r');
hold off
axis equal;

newProtocol = protocol;
newProtocol.XYZ = gridXYZcombined;
newProtocol.AFoffset = AFCombined;

%% generate new protocol file

writeFusionProtocol(protocolFile, newProtocol);

