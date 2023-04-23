clear; close all; clc

%% load positions object
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
dataDir = scriptPath;
% imdir = fullfile(scriptPath,'imgs');
% imdir = 'G:\My Drive\Research\Heemskerk lab\Data\220928_live_fixed\corrected_tracking\imgs';
imdir = fullfile(dataDir,'imgs');
load(fullfile(dataDir,'positions.mat'),'positions')
load(fullfile(dataDir,'meta.mat'),'meta')

ntime = positions(1).nTime;
nucChannel = meta.nucChannel;
fnameformat = 'stitched_MIP_p%.4d_w%.4d_t%.4d.tif';
for ii = 1:length(positions)
    positions(ii).filenameFormat = fnameformat;
end

A = positions(1).loadImage(imdir,0,1);
npos = length(positions);
imsize = size(A); imsize = imsize(1:2);

saveDir = scriptPath;

if exist(fullfile(saveDir,'validatedTracking.mat'),'file')
    load(fullfile(saveDir,'validatedTracking.mat'))
else
    UL = NaN(npos,2);
    P(npos) = Position(meta);
    verified = cell(1,npos);
    imsizes = NaN(npos,2);
end

cellfun(@numel,verified)
cellfun(@sum,verified)

%% pick a position and do automated tracking
pidx = 1;
G = P(pidx).G;

allpoints = cell2mat({P(pidx).cellData.XY}');

max_linking_distance = 45;
max_gap_distance = 70;
max_gap_time = 5;

trackopts = struct(...
    'max_linking_distance', max_linking_distance,...
    'max_gap_distance',     max_gap_distance,...
    'max_gap_time',         max_gap_time,...
    'validate_links',       false,...
    'dataDir',              imdir,...
    'dejitter',             true,...
    'imsize',               imsize);

% H = create_tracks_v6(P(pidx), trackopts);
H = create_tracks_v4(P(pidx), trackopts);

%% build histories for verified tracks

fields = {'XY','NCratio','nucLevel','cytLevel','nucArea','nucZ',...
    'nucMajorAxis','nucMinorAxis','nucCircularity'};

%ground-truth
h1 = graphSignalingHistories(P(pidx),G,fields,find(verified{pidx}));
%automated
h2 = graphSignalingHistories(P(pidx),H,fields,find(verified{pidx}));

%%
close all
tvec = (0:ntime-1)*10/60;

first_times = cellfun(@(x) x(1), {h2.Time});

idxs = find(first_times == 1);
d = zeros(length(idxs),1);
S1 = zeros(length(idxs),ntime); S2 = S1;
for ii = 1:length(idxs)
    S1(ii,h1(idxs(ii)).Time) = h1(idxs(ii)).CellIdxs;
    S2(ii,h2(idxs(ii)).Time) = h2(idxs(ii)).CellIdxs;
    d(ii) = pdist2(S1(ii,:),S2(ii,:),'jaccard');
    
%     s1 = S1(ii,:); s2 = S2(ii,:);
%     diff = double(abs(s1-s2)>0); diff(s1==0 & s2==0) = NaN;
%     plot(tvec,diff); ylim([-0.5,1.5]); xlim([min(tvec),max(tvec)])
%     cleanSubplot
%     title(sprintf('Percent incorrect = %g',100*d(ii)))
%     pause
end
% 
close all
histogram(d)
xlabel('Jaccard index for full tracks (% incorrect)')
cleanSubplot

fprintf('\n# of correct tracks = %d\n',sum(d<=0.15))
fprintf('# of incorrect tracks = %d\n',sum(d>0.15))

% linkingAccuracy(G,H,h1);
Dout = outdegree(G); DoutH = outdegree(H);
[accuracy, pred1, pred2] = linkingAccuracy(G,H,h1);
fp = sum(pred1 == 0 & pred2 > 0); %false positives
fn = sum(pred1 > 0 & pred2 == 0); %false negatives

swaps = sum(pred1 > 0 & pred2 > 0 & ~(pred1 == pred2));

badlinks = find(abs(pred1 - pred2) > 0);
badpreds = unique([pred1(badlinks),pred2(badlinks)]);
badpreds = badpreds(badpreds > 0);

added_divs = sum(DoutH(badpreds) == 2);
missed_divs = sum(Dout(badpreds) == 2);

fprintf('# false negatives = %d\n',fn)
fprintf('# false positives = %d\n',fp)
fprintf('# swaps = %d\n',swaps)
fprintf('# missed cell divisions = %d\n',missed_divs)
fprintf('# false positive cell divisions = %d\n',added_divs)

figure
histogram(tvec(G.Nodes.frame(badpreds)))
xlabel('Time (hours)')
ylabel('Frequency')
title('Distribution of time points at which tracking errors occur')
cleanSubplot
