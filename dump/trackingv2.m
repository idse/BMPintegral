clear all; close all;

%addpath(genpath('/Users/idseimac/stemcells')); 
%addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 

%[scriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);

dataDir = '/Users/idse/data_tmp/170705_Gridse_reduced';
MIPdir = fullfile(dataDir,'MIP');
meta = MetadataAndor(dataDir);

%% read tracking data

moves = {};
stats = {};
splits = {};
tmax = 10; %meta.tPerFile;
labels = zeros([1024 1024 tmax],'uint16');

for ti = 1:tmax
    
    fname = sprintf('%0.5d.h5',ti-1);
    fullfname = fullfile(MIPdir, 'Gridse_MIP_p0000_w0001_H5-Event-Sequence.h5', fname);
    info = h5info(fullfname,'/tracking');

    %dataSets = [info.Groups.Datasets];
    if ~isempty(info.Datasets)
        gotMoves = ~isempty(strmatch('Moves', {info.Datasets.Name}));
        gotSplits = ~isempty(strmatch('Splits', {info.Datasets.Name}));
    else
        gotMoves = false;
        gotSplits = false;
    end

    trackID{ti} = h5read(fullfname, '/objects/meta/id');
    if gotMoves
        moves{ti} = h5read(fullfname, '/tracking/Moves');
    end
    if gotSplits
        splits{ti} = h5read(fullfname, '/tracking/Splits');
    end
    labels(:,:,ti) = h5read(fullfname, '/segmentation/labels');
    stats{ti} = regionprops(labels(:,:,ti),'Area','Centroid');
    
    % progress indicator
    fprintf('.');
    if mod(ti,80)==0
        fprintf('\n');
    end
end
fprintf('\n');

%%

% right now extractData uses CC from nucmask which is cleaned up
% with Ilastik tracking we can skip cleanup
% but need to make PixelIdxList from tracking labelmatrix
% we can just do this with regionprops

ti = 1;
imshow(labels(:,:,ti),[])
CC = bwconncomp(labels(:,:,ti)>1);

%% chain together the moves

tracks = {};
A = {};
X = {};
trackSplits = [];

% initialize tracks
for ti = 1:2
    A{ti} = cat(1,stats{ti}.Area);
    X{ti} = cat(1,stats{ti}.Centroid);
end
% tracks.idx is a chain of object indices
for j = 1:size(moves{2},2)
    tracks{j}.idx = moves{2}(:,j)';
    tracks{j}.A = [A{1}(tracks{j}.idx(1)) A{2}(tracks{j}.idx(2))];
    tracks{j}.X = cat(1,X{1}(tracks{j}.idx(1),:), X{2}(tracks{j}.idx(2),:));
end

for ti = 3:tmax
    
    tracksStructArr = cat(1,tracks{:});
    tracksMat = cat(1,tracksStructArr.idx);
    % a vector of object indices at the end of tracks (at previous time)
    tracksEnd = tracksMat(:,end);
    
    A = cat(1,stats{ti}.Area);
    X = cat(1,stats{ti}.Centroid);
    
    for j = 1:size(moves{ti},2)
        
        % match the moves to the track ends
        tri = tracksEnd == moves{ti}(1,j);
        l = moves{ti}(2,j);
        
        % if match found and this is the first match 
        % (determined by the track length being shorter than current time)
        if any(tri) && numel(tracks{tri}.idx) < ti
            
            tracks{tri}.idx = [tracks{tri}.idx l];
            tracks{tri}.A = [tracks{tri}.A A(l)];
            tracks{tri}.X = cat(1,tracks{tri}.X, X(l,:));
            
        % if split (division, so track already has entry at this time) : 
        % duplicate track and update end to other daughter
        elseif any(tri) && numel(tracks{tri}.idx) == ti 
            
            tracks{end+1} = tracks{tri};
            tracks{end}.idx(end) = l;
            tracks{end}.A(end) = A(l);
            tracks{end}.X(end,:) = X(l,:);
        
            trackSplits = cat(1, trackSplits, [ti find(tri) numel(tracks)]);
            
        % if no match, start new track
        else
            tracks{end+1} = struct(...
                    'idx', [zeros([1 ti-1]) l],...
                    'A', A(l),...
                    'X', X(l,:));
        end
    end
    % for tracks where no match was found, add 0 idx to the end
    % to make this matrix of indices all need to be of the same lenght
    for tri = 1:numel(tracks)
       if numel(tracks{tri}.idx) < ti
           tracks{tri}.idx = [tracks{tri}.idx 0];
       end
    end 
end

tracksStructArr = cat(1,tracks{:});
tracksMat = cat(1,tracksStructArr.idx);
realTracks = tracks(sum(tracksMat > 0,2) > 1);

%% visualize time point

figure,
ti = 3;
X = cat(1, stats{ti}.Centroid);
textlabels = strsplit(num2str(1:size(X,1)),' ');

imshow(labels(:,:,ti),[])
hold on
text(X(:,1),X(:,2),textlabels,'Color','g');
hold off

%% visualize tracks

ti = 3;

%colors = lines(numel(realTracks));
colors = jet(tmax);
tracksStructArr = cat(1,tracks{:});
tracksMat = cat(1,tracksStructArr.idx);
%trackSelection = find(tracksMat(:,1)> 0);
trackSelection = 1:numel(tracks);

clf
imshow(labels(:,:,ti),[])
hold on
for tri = trackSelection
    x = realTracks{tri}.X(:, 1);
    y = realTracks{tri}.X(:, 2);
    p = line(x, y, 'LineWidth',2);
    colorsP = colors(realTracks{tri}.idx > 0,:);
    colorsP = [255*colorsP x*0+255];
    colorsP = uint8(colorsP');
    drawnow
    set(p.Edge, 'ColorBinding','interpolated', 'ColorData',colorsP)
end
ti = 133;
%scatter(realTracks{tri}.X(ti, 1),realTracks{tri}.X(ti, 2),'MarkerEdgeColor','r');%,'MarkerSize',100)
hold off

%% single cell signal histories




