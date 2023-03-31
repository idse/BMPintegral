function makeSignalingGraphs(positions, meta, options)

%     'wellnrs', 1, 
%     'channel', S4Channel,...
%     'mode','N:C',...
%     'minNCells', 10,...
%     'dataChannels', opts.dataChannels,...
%     'treatmentTime', treatmentTime,...
%     'signalingRange', [0.6 1.2]
%     'yLabel',         'ERKKTR'

    if ~isfield(options,'mode')
        options.mode = 'N:C';
    end
    SNRcutoff = ones([meta.nChannels 1]);
    if isfield(options,'SNRcutoff')
       SNRcutoff(options.channel) = options.SNRcutoff;
    end
    if ~isfield(options,'use_distinguishable_colors')
        options.use_distinguishable_colors = false;
    end
    if ~isfield(options,'excludeTimes')
        excludeTimes = [];
    else
        excludeTimes = options.excludeTimes;
    end

    wellnrs = options.wellnrs;

    channelIdx = find(options.dataChannels == options.channel);
    
    % experiment parameters
    %--------------------

    posPerCondition = meta.posPerCondition;
    
    % time
    tmax = meta.nTime;
    s = strsplit(meta.timeInterval,' ');    
    dt = str2double(s{1});
    unit = s{2};
    t = ((1:tmax) - options.treatmentTime)*dt;

    % display parameters
    %--------------------
    
    % color
    fgc = 'k';
    bgc = 'w';
    graphbgc = 1*[1 1 1]; 
    graphfgc = 'r';
    colors = hsv(posPerCondition); %0.5*[1 0 0]
    %colMultiConditions = jet(numel(wellnrs));
    if options.use_distinguishable_colors
        colorMultiConditions = distinguishable_colors(numel(wellnrs));
    else
        colorMultiConditions = lines(numel(wellnrs));
    end
    
    % axes
    xrange = [t(1), t(end)];
    yrange = options.signalingRange;

    clf
    hold on
        
    for wi = 1:numel(wellnrs)
        
        wellnr = wellnrs(wi);

        % calculate means over positions
        %-----------------------------

        if meta.loop4well
            flipwell = 9-wellnr;
            conditionPositions = [(posPerCondition/2)*(wellnr-1)+1:(posPerCondition/2)*wellnr...
                                  (posPerCondition/2)*(flipwell-1)+1:(posPerCondition/2)*flipwell];
        else
            conditionPositions = posPerCondition*(wellnr-1)+1:posPerCondition*wellnr;
        end

%         ttraceCat = cat(1,positions.timeTraces);
%         ttraceCat = ttraceCat(conditionPositions);
% 
%         % weigh by number of cells, tends to make little difference
%         W = cat(1,positions.ncells); 
%         W = W(conditionPositions,:)';
%         %W = ones(size(W)); % don't weigh
%         W = bsxfun(@rdivide, posPerCondition*W, sum(W,2));
% 
%         nucTrace = cat(3,ttraceCat.nucLevelAvg);
%         cytTrace = cat(3,ttraceCat.cytLevelAvg);
%         bgTrace = cat(3,ttraceCat.background);
% 
%         nucTrace = squeeze(nucTrace(:,channelIdx,:));
%         cytTrace = squeeze(cytTrace(:,channelIdx,:));
%         bgTrace = squeeze(bgTrace(:,channelIdx,:));

%         nucMean = nanmean(nucTrace(1:tmax,:).*W,2);
%         cytMean = nanmean(cytTrace(1:tmax,:).*W,2);
%         bgMean = nanmean(bgTrace(1:tmax,:).*W,2);
% 
%         %ratioMean = (cytMean - bgMean)./(nucMean - bgMean);
%         ratioMean = (nucMean - bgMean)./(cytMean - bgMean);

        % make graph
        %-----------------------------
        
        F = zeros([numel(conditionPositions) tmax]);
        for i = 1:numel(conditionPositions)

            pi = conditionPositions(i);

%                 nucTrace = positions(pi).timeTraces.nucLevelAvg;
%                 bgTrace = positions(pi).timeTraces.background;
%                 cytTrace = positions(pi).timeTraces.cytLevelAvg;

            %R = (cytTrace - bgTrace)./(nucTrace - bgTrace);
            %R = (nucTrace - bgTrace)./(cytTrace - bgTrace);
            
            positions(pi).makeAvgTimeTraces(SNRcutoff);
            
            if strcmp(options.mode,'N:C')
                f = positions(pi).timeTraces.ratioAvg(:,channelIdx);
            elseif strcmp(options.mode,'C:N')
                f = 1./positions(pi).timeTraces.ratioAvg(:,channelIdx);
            elseif strcmp(options.mode, 'N')
                f = positions(pi).timeTraces.nucLevelAvg(:,channelIdx);
            elseif strcmp(options.mode, 'C')
                f = positions(pi).timeTraces.cytLevelAvg(:,channelIdx);
            else
                error('unknown mode');
            end
            
            f(excludeTimes) = NaN;
            F(pi,:) = f(1:tmax)';
            F(pi, positions(pi).ncells < options.minNCells) = NaN;
        end

        if numel(wellnrs) == 1
            
            for i = 1:numel(conditionPositions)
                plot(t,F(conditionPositions(i),:),'Color', colors(i,:))
                %legend({'1','2','3','4','mean'});
                color = 'k';
            end
        else
            color = colorMultiConditions(wi,:);
        end
        FMean = nanmean(F(conditionPositions,:),1);
        plot(t, FMean, graphfgc,'LineWidth',2, 'Color', color)
    end

    fs = 24;
    xlabel(['time (' unit ')'], 'FontSize',fs,'FontWeight','Bold','Color',fgc)
    ylabel(options.yLabel, 'FontSize',fs,'FontWeight','Bold','Color',fgc);
    
    xlim(xrange);
    if ~isempty(yrange)
        ylim(yrange);
    end
    set(gcf,'color',bgc);
    set(gca, 'LineWidth', 2);
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
    set(gca,'XColor',fgc);
    set(gca,'YColor',fgc);
    set(gca,'Color',graphbgc);
    hold off 
    
end