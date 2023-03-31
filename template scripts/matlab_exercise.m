%% an image is a function
%----------------------------------------------------------------

% define a grid of x and y values
[X,Y] = meshgrid(1:1024, 1:1024);
% scale the X axis to values from 0 to 2pi
X = X*(2*pi./1024);
% define an intensity function
I = 1 + sin(2*X);
% plot that function vs x for a given y index
yi = 1;
plot(X(yi,:),I(yi,:))

% see that function as an image
figure, imshow(I,[])

%% data types and indexing
%----------------------------------------------------------------

% matrix
M = rand(3)
M(1)
M(1,2)
M(1,:)
M(:,1)
M(:)
class(M)

% math with a matrix (threshold)
M > 0.5
class(M > 0.5)

M2 = uint16(M)
class(M2)

%% cell array
C = cell(4)
C{1} = M
C(1)
C{1}

%% structure
S = struct('field1', M, 'field2', false, 'field3', {C})

%% structure array
S = struct('field1', M, 'field2', false, 'field3', C)
S = struct('field1', []);
for i = 1:3
    S(i).field1 = i;
end
S

%% basic programming
%----------------------------------------------------------------

for i = 1:10
   disp(i); 
end

threshold = 0.5;
if any(M > threshold)
    disp([num2str(sum(M(:) > threshold)) ' elements are greater than ' num2str(threshold)]);
else
    disp('all elements below threshold');
end
    
%% manipulating data

for i = 1:numel(C)
    C{i} = rand(3);
end

% concatenating
cat(2,C{1:2})
[C{1} C{2}]
[C{1} C{2}(:,1)]
    
%% reading data: define where the data is
%----------------------------------------------------------------

scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);

% raw data location, modify if not the same as the location of this script
dataDir = fullfile(scriptPath); 
cd(dataDir);

% make an array of all the filenames
listing = dir(fullfile(dataDir,'MIP','*w0001.tif'));

fi = 1;

fname = fullfile(dataDir,'MIP', listing(fi).name)
imnuc = imread(fname);
class(imnuc)
imfinfo(fname)
imshow(imnuc)

%% scale contrast display
imshow(imnuc,[])

%% scale contrast 

imdouble = mat2gray(imnuc);
tol = 0.01;
Ilim = stretchlim(imdouble, tol);
adjustedimnuc = imadjust(imdouble, Ilim);
imshow(adjustedimnuc);

%% reading the segmentation

listing = dir(fullfile(dataDir,'MIP','*h5'));
listing(1)
fname = fullfile(dataDir,'MIP',listing(fi).name);
info = h5info(fname)
info.Datasets

nuclearSegmentation = h5read(fname,'/exported_data');
size(nuclearSegmentation)
nuclearSegmentation = squeeze(nuclearSegmentation); % remove singleton dimensions
size(nuclearSegmentation)
class(nuclearSegmentation)
min(nuclearSegmentation(:))
max(nuclearSegmentation(:))

imshow(nuclearSegmentation,[])

%% logical operations and making the mask type logical

imshow(nuclearSegmentation > 1)
imshow(nuclearSegmentation == 1)
nuclearSegmentation = nuclearSegmentation == 1;
class(nuclearSegmentation)

%% overlay segmentation with image

overlay = cat(3, adjustedimnuc, adjustedimnuc + 0.3*nuclearSegmentation, adjustedimnuc);
imshow(overlay);

%% transpose: swap rows and columns

nuclearSegmentation = nuclearSegmentation';
overlay = cat(3, adjustedimnuc, adjustedimnuc + 0.3*nuclearSegmentation, adjustedimnuc);
imshow(overlay);

%% working with the nuclear mask

cleanSeg = imclearborder(nuclearSegmentation);
cleanSeg = imfill(cleanSeg,'holes');
overlay = cat(3, adjustedimnuc, adjustedimnuc + 0.3*cleanSeg, adjustedimnuc);
imshow(overlay)

%%

CC = bwconncomp(nuclearSegmentation)
% doc regionprops
stats = regionprops(CC,'Area','Solidity');
areas = cat(1, stats.Area);
hist(areas, 100)

small = areas < 200;
sum(small)

smallPixelIdx = cat(1, CC.PixelIdxList{small});
cleanerSeg = cleanSeg;
cleanerSeg(smallPixelIdx) = false;

overlay = cat(3, adjustedimnuc  + 0.3*cleanerSeg, adjustedimnuc + 0.3*cleanSeg, adjustedimnuc);
imshow(overlay)


%%

solidity = cat(1, stats.Solidity);
hist(solidity, 100)
solid = solidity > 0.9;

notsolidPixelIdx = cat(1, CC.PixelIdxList{~solid});
cleanerSeg = cleanSeg;
cleanerSeg(notsolidPixelIdx) = false;

overlay = cat(3, adjustedimnuc  + 0.3*cleanerSeg, adjustedimnuc + 0.3*cleanSeg, adjustedimnuc);
imshow(overlay)


%% combine selections with logical operation

good = solid & ~small;
badPixelIdx = cat(1,CC.PixelIdxList{~good});

cleanerSeg = cleanSeg;
cleanerSeg(badPixelIdx) = false;

overlay = cat(3, adjustedimnuc  + 0.3*cleanerSeg, adjustedimnuc + 0.3*cleanSeg, adjustedimnuc);
imshow(overlay)

%% reading out intensities within the nuclear mask

nuclearIntensities = zeros([sum(good) 1]);
goodIndex = find(good);
for i = 1:sum(good)
   nuclearIntensities(i) = mean(imnuc(CC.PixelIdxList{goodIndex(i)}));
end
hist(nuclearIntensities,50)

%% reading out brachyury reporter intensities

nucfiles = dir(fullfile(dataDir,'MIP','*MIP_*w0001.tif'));
brafiles = dir(fullfile(dataDir,'MIP','*w0000.tif'));
segfiles = dir(fullfile(dataDir,'MIP','*h5'));

fname = fullfile(dataDir,'MIP', brafiles(fi).name);
imbra = imread(fname);

im = cat(3,imbra,imnuc);

for i = 1:sum(good)
    for ci = 1:nChannels
       imc = im(:,:,ci);
       nuclearIntensities(i,ci) = mean(imc(CC.PixelIdxList{goodIndex(i)}));
    end
end

%hist(nuclearIntensities(:,1),50)
bins = 650:10:1200;
hist(nuclearIntensities(:,1),bins)
xlim([bins(1) bins(end)]);

%% comparing brachyury reporter intensities between images

% set intensity limits by brightest image
fname = fullfile(dataDir,'MIP', brafiles(19).name);
imbra = imread(fname);
Ilim = stretchlim(imbra);

allNuclearIntensities = {};
nChannels = 2;

for fi = 1:numel(brafiles)

    % read bra channel MIP
    fname = fullfile(dataDir,'MIP', brafiles(fi).name);
    imbra = imread(fname);
    imshow(imadjust(imbra, Ilim))

    % corresponding nuclear image
    fname = fullfile(dataDir,'MIP', nucfiles(fi).name);
    imnuc = imread(fname);

    % concatenate
    im = cat(3,imbra,imnuc);

    % segmentation
    fname = fullfile(dataDir,'MIP', segfiles(fi).name);
    nuclearSegmentation = h5read(fname,'/exported_data');
    nuclearSegmentation = squeeze(nuclearSegmentation==1)';
    nuclearSegmentation = imclearborder(nuclearSegmentation);
    
    CC = bwconncomp(nuclearSegmentation);
    stats = regionprops(CC,'Area','Solidity');
    areas = cat(1, stats.Area);
    solidity = cat(1, stats.Solidity);
    good = solidity > 0.9 & areas > 200;
    goodIndex = find(good);
    
    badPixelIdx = cat(1,CC.PixelIdxList{~good});
    nucSegClean = nuclearSegmentation;
    nucSegClean(badPixelIdx) = false;
    imshow(nucSegClean);

    nuclearIntensities = zeros([numel(goodIndex) nChannels]);

    for i = 1:numel(goodIndex)
        for ci = 1:nChannels
           imc = im(:,:,ci);
           nuclearIntensities(i,ci) = mean(imc(CC.PixelIdxList{goodIndex(i)}));
        end
    end
    %hist(nuclearIntensities(:,1),50)
    bins = 650:10:1200;
    hist(nuclearIntensities(:,1),bins)
    xlim([bins(1) bins(end)]);
    
    [fi numel(goodIndex)]
    
    allNuclearIntensities{fi} = nuclearIntensities;
end

%%

meanNucI = zeros([numel(nucfiles) 2]);
for fi = 1:numel(nucfiles)
    meanNucI(fi,:) = mean(allNuclearIntensities{fi});
    N = size(allNuclearIntensities{fi},1);
    errorNucI(fi,:) = std(allNuclearIntensities{fi})/sqrt(N);
end

errorbar(meanNucI(:,1),errorNucI(:,1))

%% try segment Bra in Ilastik for counting

imshow(medfilt2(imbra,[5 5]) > 900)
