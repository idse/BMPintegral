clear all; close all;

addpath(genpath('/Users/idseimac/stemcells')); 
addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 

%[scriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);

dataDir = '/Users/idse/data_tmp/170705_Gridse';
MIPdir = fullfile(dataDir,'MIP');
meta = MetadataAndor(dataDir);

%% read tracking data

moves = {};
stats = {};
splits = {};
tmax = meta.tPerFile;
labels = zeros([1024 1024 tmax],'uint16');

for ti = 1:tmax
    
    fname = sprintf('%0.5d.h5',ti-1);
    fullfname = fullfile(MIPdir, 'Gridse_MIP_p0000_w0001_H5-Event-Sequence', fname);
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
fullfname = fullfile(MIPdir, 'Tracking.h5');
bla = h5read(fullfname,'/exported_data');
    

%% chain together the moves

tracks = {};
A = {};
X = {};

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
            tracks{end}.A = A(l);
            tracks{end}.X = X(l,:);
        
        % if no match, start new track
        else
            tracks{end+1} = struct(...
                    'idx', [zeros([1 ti-1]) l],...
                    'A', [zeros([1 ti-1]) A(l)],...
                    'X', cat(1,zeros([ti-1 2]),X(l,:)));
        end
    end
    % for tracks where no match was found, add 0 to the end
    for tri = 1:numel(tracks)
       if numel(tracks{tri}.idx) < ti
           tracks{tri}.idx = [tracks{tri}.idx 0];
           %tracks{tri}.A = [tracks{tri}.A 0];
           %tracks{tri}.X = cat(1,tracks{tri}.X, [0 0]);
       end
    end 
end

tracksStructArr = cat(1,tracks{:});
tracksMat = cat(1,tracksStructArr.idx);
realTracks = tracks(sum(tracksMat > 0,2) > 1);

%% visualize time point

figure,
ti = 1;
X = cat(1, stats{ti}.Centroid);
textlabels = strsplit(num2str(1:size(X,1)),' ');

imshow(labels(:,:,ti),[])
hold on
text(X(:,1),X(:,2),textlabels,'Color','g');
hold off

%% visualize tracks

figure,
ti = 3;

colors = lines(numel(realTracks));

tracksStructArr = cat(1,tracks{:});
tracksMat = cat(1,tracksStructArr.idx);
%trackSelection = find(tracksMat(:,1)> 0);
trackSelection = 1:2:80;

imshow(labels(:,:,ti),[])
hold on
for tri = 1:numel(realTracks)
    if any(trackSelection == tracksMat(tri,1))
        x = realTracks{tri}.X(realTracks{tri}.idx > 0, 1);
        y = realTracks{tri}.X(realTracks{tri}.idx > 0, 2);
        line(x, y, 'Color',colors(tri,:),'LineWidth',2);
    end
end
hold off

% trackIdx is too long for trackPosition
