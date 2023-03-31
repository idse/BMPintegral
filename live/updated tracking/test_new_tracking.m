%Test the tracking algorithm with either real or simulated data
clear; clc;

%% Set parameters
%Parameters, constants
%max frame-frame cell movement distance in pixels; farther cells are ignored
max_linking_distance = 30; %~25-30 for 40x data; ~15 for 20x data
max_gap_distance = 70; %approximately double the max linking distance
max_gap_time = 5;

trackopts = struct(...
    'max_linking_distance', 30,...
    'max_gap_distance',     70,...
    'max_gap_time',         5,...
    'validate_links',       true);

%% Select data
%You are prompted to pick a positions.mat file
cd('Y:\Seth') %this is just to help find my data faster
[selectfile, selectpath] = uigetfile('*.mat','Select positions file');
load(fullfile(strcat(selectpath, selectfile)));
P = positions(1);
trackopts.dataDir = selectpath;

try
    img = P.loadImage(selectpath, P.nucChannel, 1);
    imsize = [size(img,1), size(img,2)];
    trackopts.imsize = imsize;
catch
    try
        img = P.loadImage(fullfile(selectpath,'MIP'), P.nucChannel, 1);
        imsize = [size(img,1), size(img,2)];
        trackopts.imsize = imsize;
        trackopts.dataDir = fullfile(selectpath,'MIP');
    catch
        warning('Could not load image: using image size of [1024, 1024]')
        imsize = [1024, 1024];
    end
end

%% Perform tracking with old(er) function
P.makeTracks(max_linking_distance,...
    max_gap_distance, max_gap_time, imsize);


%% make a new field in cellData and identify (potentially) dividing cells

for ti = 1:P.nTime
    majors = P.cellData(ti).nucMajorAxis;
    minors = P.cellData(ti).nucMinorAxis;
    area = P.cellData(ti).nucArea;
    axisRatios = majors./minors;
    divs = axisRatios >= 2.5;
    divs = divs.*(area < 750);
    P.cellData(ti).divs = divs;
end

%% Create tracks with newer tracking code
c1 = clock;
trees = create_tracks_v2(P, trackopts);
P.tree = trees;
c2 = clock;
time = c2 - c1;
time = (time(6) + time(5)*60 + time(4)*3600);
disp(strcat("Elapsed time is ", num2str(time)," seconds"))

% positions = P;
% save(fullfile(selectpath,selectfile), 'positions');

%% Visualize tracks
close all

visopts = struct('tail_length', 25,...
    'write_to_tif', true,...
    'write_to_avi', false,...
    'use_tif', true,...
    'use_seg', false,...
    'default_seg', true,...
    'default_name', true,...
    'visible', 'off');
%Making track visualization
tic
visualizeTracks4(P, visopts);
toc

%% Making trace visualization
traceopts = struct('tail_length', 25,...
    'write_to_tif', false,...
    'write_to_avi', true,...
    'use_tif', true,...
    'use_seg', false,...
    'default_seg', true,...
    'default_name', true,...
    'visible', 'off');
minlength = 90;
P.makeTimeTraces(minlength);
ntraces = length(P.timeTraces.trackT); %do all traces
channel = 0;
tic
visualizeCellSignaling(P, ntraces, channel, traceopts);
toc

%% make visualization of longest traces that make it to the end
traceopts = struct('tail_length', 25,...
    'write_to_tif', true,...
    'write_to_avi', false,...
    'use_tif', true,...
    'use_seg', false,...
    'default_seg', true,...
    'default_name', true,...
    'visible', 'off');
minlength = 10;
P.makeTimeTraces(minlength);
ntraces = 10;
channel = 0;
tic
visualizeBestTracks(P, ntraces, channel, traceopts);
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
