function [P,xi,Pstd] = radialPositive(stats, positions, condi, meta, combos,order)
    % radial probability profile
    %
    % [P,xi] = radialProbability(stats, positions, condi, meta, combos)
    % 
    % condi : index of condition/well
    % combos = {[2 3],..} for expression two or more markers (more not
    % implemented)

    bw = 15; % bandwidth for kde

    if ~exist('order','var')
        order = 1:numel(stats.markerChannels);
    end
    channels = stats.markerChannels(order);
    
    colRadius = positions((condi-1)*meta.posPerCondition + 1).radiusMicron;
    
    % part of true radius estimate below
    celldensity = {};
    R = sqrt(sum(stats.XY{condi}.^2,2))*meta.xres;
    nominalradius = stats.radiusMicron{condi};
    %-------
    
    Psample = [];
    sampleid = unique(stats.sample{condi});
    Nsamples = numel(sampleid);
	ii = 0; % keep track of total iterations to reuse xi from first
    
    for si = 1:Nsamples+1 % last iteration combines all samples
        
        if si <= Nsamples
            sampleidx = stats.sample{condi} == sampleid(si);
        else
           sampleidx = true(size(stats.sample{condi}));
        end
        Rs = R(sampleidx);
        
        %------------------------------------------------------------------------
        % estimate true radius
        %------------------------------------------------------------------------
        if si <= Nsamples
        
            % estimate cell density and colony radius (which may differ
            % slightly from nominal radius of micropattern)
            PI = 3.1415;

            bw = 5;
            [f,xip] = ksdensity(Rs,'Bandwidth',bw,'BoundaryCorrection','reflection');

            dx = xip(2)-xip(1);
            Ncells = numel(Rs);
            celldensity{si} = Ncells*f./(2*PI*xip);

            xrangeinside = xip > bw & xip < 0.9*nominalradius;
            avgdensityinside = sum(2*PI*celldensity{si}(xrangeinside).*xip(xrangeinside)*dx)/sum(2*PI*xip(xrangeinside)*dx);
            % define the radius as the point where the density drops to 10% of
            % the mean inside
            fullxrange = xip > bw & celldensity{si} > 0.1*avgdensityinside;
            trueradius = max(xip(fullxrange));

            %if options.useTrueRadius
                radius = trueradius;
                disp('adjusting for true radius of colony');
            %else
    %            radius = nominalradius;
    %         end

            %convert radius to edge distance
            Rs = radius - Rs;
            R(sampleidx) = Rs;
        end
   
        for ci = channels
            
            posidx = stats.nucLevel{condi}(sampleidx,ci) > stats.thresholds(ci);
            posR = Rs(posidx);

            % histogram approach
            %[N,edges] = histcounts(R);
            %[Npos,~] = histcounts(posR,edges);
            %xi = edges(1:end-1) + diff(edges)/2;
            %bar(xi, Npos./N)
            %plot(xi, Npos./N)
            %nansum(Npos./N)

            xi = linspace(0,nominalradius,100);
            if ii == 0
                P = zeros(max(channels), numel(xi));
            end
            [f,~] = ksdensity(Rs,xi,'BoundaryCorrection','reflection','Bandwidth',bw);%,'Support','positive');
            %dx = diff(xi);
            %dx = dx(1);
            %intf = dx*sum(f);
            if sum(posidx) > 0
                [fpos,~] = ksdensity(posR,xi,'BoundaryCorrection','reflection','Bandwidth',bw);%,'Support','positive');
            else
                fpos = xi*0;
                warning(['no positive cells for ' meta.channelLabel{ci}]);
            end
            %intfpos = dx*sum(fpos);

            % normalize
            P(ci,:) = (sum(posidx)/sum(sampleidx))*fpos./f;
            
            ii = ii + 1;
        end
        
        if si == 1
            Psample = P;
        else
            Psample = cat(3,Psample,P);
        end
    end

    Pstd = zeros(size(Psample(:,:,1)));
    Pavg = zeros(size(Psample(:,:,1)));
    for ci = 1:4
        Pstd(ci,:) = std(Psample(ci,:,1:Nsamples),[],3);
        Pavg(ci,:) = mean(Psample(ci,:,1:Nsamples),3);
    end

    if exist('combos','var')
        
        Pc = zeros(numel(combos), numel(xi));

        Pcsample = [];
        sii = unique(stats.sample{condi});
        
        for si = 1:Nsamples+1 % last iteration combines all samples
            
            if si == Nsamples+1
               sampleidx = true(size(stats.sample{condi}));
            else
                sampleidx = stats.sample{condi} == sii(si);
            end
            
            for ci = 1:numel(combos)

                posidx = stats.nucLevel{condi}(:,combos{ci}(1)) > stats.thresholds(combos{ci}(1)) & ...
                            stats.nucLevel{condi}(:,combos{ci}(2)) > stats.thresholds(combos{ci}(2));
                posidx = posidx & sampleidx;
                
                if sum(posidx) > 0
                    XY = stats.XY{condi};
                    posXY = XY(posidx,:);

                    Rs = sqrt(XY(sampleidx,1).^2 + XY(sampleidx,2).^2);
                    posR = sqrt(XY(posidx,1).^2 + XY(posidx,2).^2);

                    [f,~] = ksdensity(Rs,xi,'BoundaryCorrection','reflection','Bandwidth',bw);%,'Support','positive');
                    [fpos,~] = ksdensity(posR,xi,'BoundaryCorrection','reflection','Bandwidth',bw);%,'Support','positive');

                    Pc(ci,:) = (sum(posidx)/sum(sampleidx))*fpos./f;
                else
                    warning(['no positive cells for combo ' num2str(combos{ci})]);
                end
            end
            
            if si == 1
                Pcsample = Pc;
            else
                Pcsample = cat(3,Pcsample,Pc);
            end
        end
        
        Pcstd = zeros(size(Pcsample(:,:,1)));
        Pcavg = zeros(size(Pcsample(:,:,1)));
        for ci = 1:numel(combos)
            Pcstd(ci,:) = std(Pcsample(ci,:,1:Nsamples),[],3);
            Pcavg(ci,:) = mean(Pcsample(ci,:,1:Nsamples),3);
        end
    end
    
    colors = lines(8);
    colors = colors([2 5 1 3 4 6 7 8],:);

    lw = 3;
    figure,
    hold on
    for i = 1:numel(channels)
        %plot(xi*meta.xres,P(channels(i),:),'LineWidth',lw,'Color',colors(i,:));
        plot(xi,P(channels(i),:),'LineWidth',lw,'Color',colors(i,:));
    end
    if exist('combos','var')
        for i = 1:numel(combos)
            %plot(xi*meta.xres,Pc(i,:),'LineWidth',lw,'Color',colors(numel(channels) + i,:));
            plot(xi,Pcavg(i,:),'LineWidth',lw,'Color',colors(numel(channels) + i,:));
        end
    end
    % error bars
    for i = 1:numel(channels)
        %errorbar(xi*meta.xres,P(channels(i),:),Pstd(channels(i),:),'LineWidth',lw,'Color',colors(i,:));
        good = ~isnan(P(channels(i),:));
        fill([xi,fliplr(xi)],...
            [P(channels(i),good) + Pstd(channels(i),good), fliplr(P(channels(i),good) - Pstd(channels(i),good))],...
            colors(i,:),'FaceAlpha',0.2,'EdgeColor','none');
    end
    if exist('combos','var')
        for i = 1:numel(combos)
            good = ~isnan(Pc(i,:));
            fill([colRadius - meta.xres*xi(good),fliplr(colRadius-meta.xres*xi(good))],[Pc(i,good) + Pcstd(i,good), fliplr(Pc(i,good) - Pcstd(i,good))],colors(i,:),'FaceAlpha',0.2,'EdgeColor','none');
        end
    end
    hold off
    
    legendstr = meta.channelLabel(channels);
    for ci = 1:numel(combos)
        x = meta.channelLabel(combos{ci});
        legendstr = [legendstr {[x{1} ' & ' x{2}]}];
    end
    fs = 20;
    legend(legendstr,'FontSize',10,'Location','NorthEast')
    axis square
    fgc = 'k';
    bgc = 'w';
    graphbgc = 1*[1 1 1]; 

    xlim([0 colRadius]);
    xlabel('edge distance ( um )', 'FontSize',fs,'FontWeight','Bold','Color',fgc)
    ylabel('positive fraction ', 'FontSize',fs,'FontWeight','Bold','Color',fgc);

    set(gcf,'color',bgc);
    set(gca, 'LineWidth', 2);
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
    set(gca,'XColor',fgc);
    set(gca,'YColor',fgc);
    set(gca,'Color',graphbgc);
end