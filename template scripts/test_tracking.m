%Test the tracking algorithm with either real or simulated data
clear; clc;

%% Set parameters
%Parameters, constants
%max frame-frame cell movement distance in pixels; farther cells are ignored
max_linking_distance = 30; %~25-30 for 40x data; ~15 for 20x data
max_gap_distance = 70; %approximately double the max linking distance
max_gap_time = 5;

%% Select data
%You are prompted to pick a positions.mat file
cd('Y:\Seth') %this is just to help find my data faster
[selectfile, selectpath] = uigetfile('*.mat','Select positions file');
load(fullfile(strcat(selectpath, selectfile)));
P = positions(1);

try
    img = P.loadImage(selectpath, P.nucChannel, 1);
    imsize = [size(img,1), size(img,2)];
catch
    warning('Could not load image: using image size of [1024, 1024]')
    imsize = [1024, 1024];
end

%% Create tracks
P1.makeTracks(max_linking_distance,...
    max_gap_distance, max_gap_time, imsize);
% positions = P;
%save the computed tracks
% save(fullfile(selectpath,'positions.mat'), 'positions')

%% Visualize tracks
close all

opts = struct('tail_length', 25,...
    'write_to_tif', false,...
    'write_to_avi', true,...
    'use_tif', true,...
    'use_seg', false,...
    'default_seg', true,...
    'default_name', true,...
    'visible', 'off');

tic
visualizeTracks4(P, opts);
toc

%% Visualize select tracks and time traces
close all
minlength = 10;
P.makeTimeTraces(minlength);
opts = struct('tail_length', 25,...
    'write_to_tif', false,...
    'write_to_avi', false,...
    'use_tif', true,...
    'use_seg', false,...
    'default_seg', true,...
    'default_name', true,...
    'visible', 'on');

visualizeCellSignaling(P, 3, opts);

%% Dejitter video
P1 = positions(1);
for ti = 1:P1.nTime
    if isempty(P1.cellData(ti).XY)
        P1.cellData(ti) = P1.cellData(ti-1);
    end
end

[selectfile, selectpath] = uigetfile('*.tif','Select image file');
fname = fullfile(selectpath,selectfile);
imf = imfinfo(fname);
ntime = length(imf);
imprev = imread(fname,1);
disps = zeros(ntime,2);
for ti = 2:ntime
    im = imread(fname,ti);
    if sum(im,'all') > 7e8
        [shiftx,shifty,~,~] = xcorr2fft(imprev, im);
        disps(ti,:) = [shiftx, shifty];
        imprev = im;
    end
end

cumdisps = cumsum(disps,1);

for ti = 1:ntime
    P1.cellData(ti).XY = P1.cellData(ti).XY - cumdisps(ti,:);
end

P1.nucChannel = 0;

%% Create tracks with newmergelinker
c1 = clock;
trees = create_tracks_v2(P1, max_linking_distance,...
max_gap_distance, max_gap_time, imsize);
P1.tree = trees;
c2 = clock;
time = c2 - c1;
time = (time(6) + time(5)*60 + time(4)*3600);
disp(strcat("Elapsed time is ", num2str(time)))

% Visualize tracks
close all

trackopts = struct('tail_length', 25,...
    'write_to_tif', false,...
    'write_to_avi', true,...
    'use_tif', true,...
    'use_seg', false,...
    'default_seg', true,...
    'default_name', true,...
    'visible', 'off');
%Making track visualization
tic
visualizeTracks4(P1, trackopts);
toc

%Making trace visualization
traceopts = struct('tail_length', 25,...
    'write_to_tif', false,...
    'write_to_avi', true,...
    'use_tif', true,...
    'use_seg', false,...
    'default_seg', true,...
    'default_name', true,...
    'visible', 'off');
minlength = 10;
P1.makeTimeTraces(minlength);
tic
visualizeCellSignaling(P, 10, traceopts);
toc

%%
positions = P1;
save(fullfile(selectpath,selectfile),'positions')

%% make visualization of longest traces that make it to the end
traceopts = struct('tail_length', 25,...
    'write_to_tif', false,...
    'write_to_avi', true,...
    'use_tif', true,...
    'use_seg', false,...
    'default_seg', true,...
    'default_name', true,...
    'visible', 'off');
minlength = 10;
P1.makeTimeTraces(minlength);
tic
visualizeBestTracks(P1, 10, traceopts);
toc

%% junk
% Select data
% if true && use_sim_data
%     %If simulating data, set simulation options
%     %Sets roughly the percentage of cells that randomly fail to be identified
%     segment_errors = 2;
%     %Number of time points created
%     ntimepoints = 150;
%     %Controls how much cells move between frames and how much intensity
%     %changes
%     velocity = 10;
%     intensity_drift = 1;
%     %Number of cells simulated
%     numcells = 10;
%     %Cell generation area around the image window - cells are given initial
%     %locations within x=[-buffer/2,imwidth+buffer/2], y=[-buffer/2,imheight+buffer/2]
%     buffer = -300;
%     %Division frequency determines the likelihood of a cell division (this
%     %is given as a percentage)
%     division_freq = 1;
%     
%     %Generate the data
%     [P, sim_trees, sim_points] = simulateCells(numcells, ntimepoints,...
%         division_freq, segment_errors, velocity, intensity_drift,...
%         imsize, buffer);
% else
%     %If using real data, you are prompted to pick a positions.mat file
%     [selectfile, selectpath] = uigetfile('*.mat','Select positions file');
%     load(fullfile(strcat(selectpath, selectfile)));
%     P = positions(1);
%     P.nucChannel = nucChannel - 1;
% end


%compute tracks
%     trees = create_tracks(P, max_linking_distance, max_gap_distance,...
%         max_gap_time, imsize, nucChannel);
