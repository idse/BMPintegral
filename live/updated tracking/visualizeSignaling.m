function visualizeSignaling(lt, pidx, opts)
%visualize time traces for a LineageTrace object
%the 

field = opts.field;
channel = opts.channel;

P = lt.live_position(pidx);
nucChannel = P.nucChannel;
trackT = {lt.signalingHistories{pidx}.trackT};
XY = {lt.signalingHistories{pidx}.XY};
signals = {lt.signalingHistories{pidx}.(field)};
signals = cellfun(@(x) x(:,channel+1), signals, 'UniformOutput', false);
S = NaN(P.nTime,length(signals));
for ii = 1:length(signals)
    S(trackT{ii},ii) = signals{ii};
end

mintime = cellfun(@min, trackT);
mintime = min(mintime);
maxtime = cellfun(@max, trackT);
maxtime = max(maxtime);

minsignal = min(cellfun(@min, signals));
maxsignal = max(cellfun(@max, signals));

%Define an array of colors for visualization
%Distuinguishable colors function makes at most 9000 different colors; if
%there are > 9000 tracks, start repeating colors
n = numel(signals);
if n <= 9000
    colors = distinguishable_colors(n,{'g','k'});
else
    colors = distinguishable_colors(9000);
    numreps = ceil(n/9000) - 1;
    remainder = mod(n,9000);
    colors = vertcat(repmat(colors,numreps), colors(1:remainder,:));
end
%set name for saving results as a video
writepath = fullfile(opts.dataDir,'Tracks');
if ~exist(writepath,'dir')
        mkdir(writepath);
end
aviname = fullfile(writepath,[P.filename,...
    '000',num2str(nucChannel-1),'_signaling_', date, '.avi']);
%if writing as an avi video, create a VideoWriter object
if opts.write_to_avi
    disp(strcat("Writing results as video to ", aviname))
    v = VideoWriter(aviname);
    open(v)
end

f1 = figure('Units', 'normalized', 'Outerposition', [0, 0.04, 1, 0.96],...
    'Visible', opts.visible);
for ti = mintime:maxtime
    clf
    points = P.cellData(ti).XY;
    img = max(loadImage(P, opts.dataDir, nucChannel, ti),[],3);
    im2 = max(loadImage(P, opts.dataDir, channel, ti),[],3);
    if ti == mintime
        lim1 = stretchlim(img);
        lim2 = stretchlim(im2);
    end
    img = imadjust(img, lim1);
    im2 = imadjust(im2, lim2);
    
    ax1 = subplot(1,3,1);
    hold(ax1,'on')
    imshow(img)
    scatter(points(:,1),points(:,2),15,'g','filled')
    ax3 = subplot(1,3,2);
    hold(ax3,'on')
    imshow(im2)
    scatter(points(:,1),points(:,2),15,'g','filled')
    ax2 = subplot(1,3,3);
    hold(ax2,'on')
    for ii = 1:length(XY)
        startT = max(mintime, ti - opts.tail_length);
        idxs = find(trackT{ii} >= startT & trackT{ii} <= ti);
        x = XY{ii}(idxs,1);
        y = XY{ii}(idxs,2);
        if ~isempty(x)
            line(ax1,x,y,'Color',colors(ii,:),'LineWidth', 1)
            line(ax3,x,y,'Color',colors(ii,:), 'LineWidth', 1)
        end
        plot(ax2, mintime:ti, S(mintime:ti,ii), 'Color', colors(ii,:))
    end
    xlim(ax2, [mintime, maxtime])
    ylim(ax2, [minsignal maxsignal])
    xlabel(ax2,'Time (frames)')
    ylabel(ax2, field)
    title(ax2, 'Signaling')
    sgtitle(sprintf('Time = %d', ti))
    
    ax1.Position([1,3]) = [0.01, 0.32];
    ax3.Position([1,3]) = [0.34, 0.32];
    ax2.Position([1,3]) = [0.69, 0.3];
    
    if opts.write_to_avi
        frame = getframe(f1);
        writeVideo(v, frame.cdata);
        if ti == maxtime
            %Close the VideoWriter after the last frame
            close(v);
        end
    end
    
    if strcmp(opts.visible,'on')
        pause(0.05)
    end
    
end

end