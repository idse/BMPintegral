clear; close all

%% load data
dataDir = 'Y:\Seth\200323_subsampled_for_tracking';
load(fullfile(dataDir,'positions.mat'))
load(fullfile(dataDir,'meta.mat'))

ntime = positions(1).nTime;
nucChannel = positions(1).nucChannel;
tail_length = 12;

A = positions(1).loadImage(dataDir,0,1);
imsize = size(A);
clear A
%position index
pidx = 1;
%XY positions of all points in all frames in one array, where row index
%corresponds to node index in the tracking digraph
allpoints = cell2mat({positions(pidx).cellData.XY}');

%% Do single-cell tracking
max_linking_distance = 60;
max_gap_distance = 70;
max_gap_time = 3;

trackopts = struct(...
    'max_linking_distance', max_linking_distance,...
    'max_gap_distance',     max_gap_distance,...
    'max_gap_time',         max_gap_time,...
    'validate_links',       false,...
    'dataDir',              dataDir,...
    'dejitter',             false,...
    'imsize',               imsize);

c1 = clock;
%do tracking, store links between cells in a directed graph (digraph)
G = create_tracks_v5(positions(pidx), trackopts);

disp('Building tracks')
%build tracks from digraph
tracks = graphtracks(G, ntime);
%get X,Y information for tracks for visualization
[X,Y] = trackXY(tracks, allpoints);
c2 = clock;
time = c2 - c1;
time = time(end) + 60*time(end-1) + 3600*time(end-2);
disp(strcat("Total time is ", num2str(time), " seconds"))

%% load images
%load all images and adjust contrast ahead of time to speed up
%interactivity later
imgs = zeros(imsize(1), imsize(2), ntime, 'uint8');
for ti = 1:ntime
    fprintf('.')
    if mod(ti,45) == 0
        fprintf('\n')
    end
    img = positions(pidx).loadImage(dataDir,nucChannel,ti);
    img = max(img,[],3);
    if ti == 1
        %make contrast adjustment limits consistent between all frames
        lim = stitchedlim(img);
    end
    %adjust contrast first, then convert to 8-bit to save RAM
    imgs(:,:,ti) = im2uint8(imadjust(img,lim));
end
fprintf('\n')

%% Do visualization/correction
ntracks = length(tracks);
%used for converting from [time, cell] indices to graph node indices
cumtimes = [0, cumsum(positions(pidx).ncells)];

%define colors for each track
colors = make_colors(ntracks);

close all
f = figure('WindowState','maximized');
%initialize a plot object for image display
p1 = imshow(imgs(:,:,1),'InitialMagnification','fit');
hold on
%initialize a scatterplot of cell centroid positions
xys = positions(pidx).cellData(1).XY;
f1 = scatter(xys(:,1),xys(:,2),25,'g','filled');

%initialize an array of Line objects to display tracking tails
p = initialize_lines(tracks, colors);

drawnow
tidx = 1;
cleanSubplot
breakvar = false;
updateflag = false;
while breakvar == false
    %update image for new frame
    set(p1,'CData',imgs(:,:,tidx));
    %update scatter plot of cell centroid positions on the image
    xys = positions(pidx).cellData(tidx).XY;
    set(f1,'XData',xys(:,1),'YData',xys(:,2));
    
    idx = 1;
    %draw tracking tails
    for ii = 1:ntracks
        for k = 1:size(tracks{ii},2)
            sidx = max(1, tidx - tail_length);
            x = X{ii}(sidx:tidx,k);
            x = x(~isnan(x));
            y = Y{ii}(sidx:tidx,k);
            y = y(~isnan(y));
            if ~isempty(x)
                set(p(idx),'XData',x,'YData',y)
            else
                set(p(idx),'XData',1,'Ydata',1)
            end
            idx = idx + 1;
        end
    end
    
    title(sprintf('%d of %d', tidx, ntime))
    
    waitforbuttonpress;
    key = f.CurrentCharacter;
    
    if strcmp(key,'a')
        %go back one time point
        tidx = tidx - 1;
    elseif strcmp(key,'d')
        %advance one time point
        tidx = tidx + 1;
    elseif strcmp(key,'w')
        %advance ten time points
        tidx = tidx + 10;
    elseif strcmp(key,'s')
        %go back ten time points
        tidx = tidx - 10;
    elseif strcmp(key,'+')
        %zoom in 2X
        xl = get(gca,'xlim');
        yl = get(gca,'ylim');
        set(gca,'xlim', mean(xl) + 0.25*(xl(2) - xl(1))*[-1 1])
        set(gca,'ylim', mean(yl) + 0.25*(yl(2) - yl(1))*[-1 1])
    elseif strcmp(key,'-')
        %zoom out 2X
        xl = get(gca,'xlim');
        yl = get(gca,'ylim');
        set(gca,'xlim', mean(xl) + (xl(2) - xl(1))*[-1 1])
        set(gca,'ylim', mean(yl) + (yl(2) - yl(1))*[-1 1])
    elseif strcmp(key,'2')
        %pan down
        yl = get(gca,'ylim');
        set(gca,'ylim', yl + 0.25*(yl(2) - yl(1)))
    elseif strcmp(key,'8')
        %pan up
        yl = get(gca,'ylim');
        set(gca,'ylim', yl - 0.25*(yl(2) - yl(1)))
    elseif strcmp(key,'4')
        %pan left
        xl = get(gca,'xlim');
        set(gca,'xlim', xl - 0.25*(xl(2) - xl(1)))
    elseif strcmp(key,'6')
        %pan right
        xl = get(gca,'xlim');
        set(gca,'xlim', xl + 0.25*(xl(2) - xl(1)))
    elseif strcmp(key,'5')
        %reset field of view to default (full image)
        set(gca,'xlim', [0,size(imgs,1)] + [-0.5 0.5])
        set(gca,'ylim', [0,size(imgs,2)] + [-0.5 0.5])
    elseif strcmp(key,'e')
        %exit the while loop
        breakvar = true;
    elseif strcmp(key,'m')
        clear s t
        %make a change to the tracking results; add or remove a link
        [lindx, ~] = listdlg('ListString',...
            {'Add link','Remove link before','Remove link after','Cancel'});
        if lindx == 1 %add link
            %keep track of source node to add 
            s = pickpoint(xys,tidx,cumtimes,'source'); %source node
        elseif lindx == 2 %remove link before (into) selected node
            t = pickpoint(xys,tidx,cumtimes,'target'); %target node
            if indegree(G,t) == 1
                s = predecessors(G,t); %source node
                G = rmedge(G,s,t); %remove edge between these nodes
                updateflag = true;
            else
                disp('There are no incoming edges to this node')
            end
            clear s t
        elseif lindx == 3 %remove link after (out of) selected node
            s = pickpoint(xys,tidx,cumtimes,'source'); %source node
            if outdegree(G,s) == 1
                t = successors(G,s); %target node
                G = rmedge(G,s,t); %remove edge between these nodes
                updateflag = true;
            else
                disp('Node has multiple edges out, pick one of these to remove a link to instead')
            end
            clear s t
        end
    elseif strcmp(key,'n')
        %pick second point for link addition (always pick source then
        %target)
        if exist('s','var') == 1
            %pick second point for added tracking link
            t = pickpoint(xys,tidx,cumtimes,'target'); %target node
            %check if (1) the second point comes after the first in time; (2)
            %if the first point has more than one link out already, (3) if the
            %second point has more than 0 links in already
            if G.Nodes.frame(t) > G.Nodes.frame(s) && outdegree(G,s) < 2 && indegree(G,t) == 0
                EdgeTable = table([s t],1,1,'VariableNames',...
                    {'EndNodes','Weight','Cost'});
                G = addedge(G,EdgeTable);
                updateflag = true;
            else
                disp('Could not add link')
            end
        else
            disp('There is no source node chosen')
        end
    elseif strcmp(key,'c')
        %cancel link addition
        clear s
    end
    
    if updateflag
        %if tracking digraph was modified, rebuild tracks and reinitialize
        %array of lines
        delete(p)
        xlabel('Rebuilding tracks')
        drawnow
        tracks = graphtracks(G, ntime);
        ntracks = length(tracks);
        
        [X,Y] = trackXY(tracks,allpoints);
        colors = make_colors(ntracks);
        p = initialize_lines(tracks, colors);
        xlabel('')
        
        updateflag = false; %reset update flag
    end
    
    if tidx > ntime
        tidx = ntime;
    elseif tidx < 1
        tidx = 1;
    end
    
    
    
    drawnow limitrate
end
close(f)

%% save most recent updated version of the tracking digraph
save(fullfile(dataDir,['G_',date,'.mat']),'G')

%% local functions
function colors = make_colors(ntracks)
%Define an array of colors for visualization
%Distuinguishable colors function makes at most 9000 different colors; if
%there are > 9000 tracks, start repeating colors
    if ntracks <= 9000
        colors = distinguishable_colors(ntracks,{'w','k'});
    else
        colors = distinguishable_colors(9000);
        numreps = ceil(ntracks/9000) - 1;
        remainder = mod(ntracks,9000); % = ntracks - numreps*9000
        colors = vertcat(repmat(colors,numreps), colors(1:remainder,:));
    end
end

function [X,Y] = trackXY(tracks,allpoints)
    X = cell(size(tracks));
    Y = cell(size(tracks));
    for ti = 1:length(tracks)
        X{ti} = NaN(size(tracks{ti}));
        Y{ti} = NaN(size(tracks{ti}));
        nonans = ~isnan(tracks{ti});
        X{ti}(nonans) = allpoints(tracks{ti}(nonans),1);
        Y{ti}(nonans) = allpoints(tracks{ti}(nonans),2);
    end
end

function p = initialize_lines(tracks, colors)
%initialize an array of Line objects to display tracking tails
    ntotal = sum(cellfun(@(x) size(x,2),tracks));
    ntracks = length(tracks);
    p = gobjects(ntotal,1);
    idx = 1;
    for ii = 1:ntracks
        for k = 1:size(tracks{ii},2)
            p(idx) = line(gca,1,1,'Color',colors(ii,:),'LineWidth',1.5);
            idx = idx + 1;
        end
    end
end

function allidx = pickpoint(xys,tidx,cumtimes,pointtype)

xlabel(strcat("Pick ", pointtype, " point"))
h = drawpoint;
xy = h.Position; %position of drawn point
delete(h)
%find the nearest cell centroid to the drawn point
D = pdist2(xys, xy(1,:));
[~, cellidx] = min(D);
%find the overall index of this point in the tracking digraph
allidx = cumtimes(tidx) + cellidx;
xlabel('')

end








