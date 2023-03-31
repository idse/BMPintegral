function res = plotRadialProfiles_old(positions, meta, options)

    if ~isfield(options,'normalize')
        options.normalize = false;
    end
    if ~isfield(options,'nucChannels') && ~isfield(options,'cytChannels')
        error('please specify channels');
    end
    if ~isfield(options,'nucChannels')
        options.nucChannels = [];
    end
    if ~isfield(options,'cytChannels')
        options.cytChannels = [];
    end
    if ~isfield(options,'legend')
        options.legend = true;
    end
    if ~isfield(options,'colors')
        options.colors = lines(5);
        options.colors = options.colors([2 5 1 3 4],:);
    end
    if ~isfield(options,'std') %|| options.normalize == true
        options.std = false;
    end
    if ~isfield(options,'individualColonies')
        options.individualColonies = false;
    end
    if ~isfield(options,'Ilimmargin')
        options.Ilimmargin = 3;
    end
    if isfield(options,'FontSize') 
        fs = options.FontSize;
    else
        fs = 20;
    end
    
    if numel(positions) == 1
        nuc_profile = positions(1).radialProfile.NucAvgSeg;
        cyt_profile = positions(1).radialProfile.CytAvgSeg;
        options.std = false;
    else
        % average if given an array of positions
        nuctmp = zeros([numel(positions) size(positions(1).radialProfile.NucAvgSeg)]);
        cyttmp = nuctmp;
        for pi = 1:numel(positions)
            nuctmp(pi,:,:) = positions(pi).radialProfile.NucAvgSeg;
            cyttmp(pi,:,:) = positions(pi).radialProfile.CytAvgSeg;
        end

        nuc_profile = squeeze(nanmean(nuctmp,1));
        cyt_profile = squeeze(nanmean(cyttmp,1));
        
        nuc_profile_std = squeeze(nanstd(nuctmp,1)); 
        cyt_profile_std = squeeze(nanstd(cyttmp,1));
    end
    
    % the edges are very noisy, because the number of cells drops to zero,
    % exclude from determining limits
    binmarg = options.Ilimmargin;
    if isfield(options,'nuclimits')
        nuclimits = options.nuclimits;
    else
        nuclimits = [min(nuc_profile(1:end-binmarg,:))' max(nuc_profile(1:end-binmarg,:))'];
    end
    if isfield(options,'cytlimits')
        cytlimits = options.cytlimits;
    else
        cytlimits = [min(cyt_profile)' max(cyt_profile)'];
    end
    
    if options.normalize
        for ci = 1:size(nuc_profile,2)
            nuc_profile(:,ci) = (nuc_profile(:,ci) - nuclimits(ci,1))./(nuclimits(ci,2) - nuclimits(ci,1));
            cyt_profile(:,ci) = (cyt_profile(:,ci) - cytlimits(ci,1))./(cytlimits(ci,2) - cytlimits(ci,1));
            nuc_profile_std(:,ci) = nuc_profile_std(:,ci)./(nuclimits(ci,2) - nuclimits(ci,1));
            cyt_profile_std(:,ci) = cyt_profile_std(:,ci)./(cytlimits(ci,2) - cytlimits(ci,1));
        end
    end
    r = positions(1).radialProfile.BinEdges(1:end-1)*meta.xres;
    %r = (r(1:end-1)+r(2:end))/2;
    
    r = r(end) - r;
    %fpositions(1).radiusMicron - r;

    %figure,
    
    legendentries = {};
    binmarg = 0;
%     clf
    hold on
    if ~isempty(options.nucChannels)

      
        for cii = numel(options.nucChannels):-1:1
            ci = options.nucChannels(cii) + 1;
            
            plot(r(1:end-binmarg), nuc_profile(1:end-binmarg,ci),'LineWidth',3, 'Color', options.colors(cii,:))

            if isempty(options.cytChannels)
                prefix = [];
            else
                prefix = 'nuc ';
            end
            legendentries = [legendentries, [prefix meta.channelLabel{ci}]];
        end
        for cii = numel(options.nucChannels):-1:1
            ci = options.nucChannels(cii) + 1;
            if options.std
                errorbar(r(1:end-binmarg), nuc_profile(1:end-binmarg,ci), nuc_profile_std(1:end-binmarg,ci),'LineWidth',1, 'Color', options.colors(cii,:)); 
            end
        end
        
        if options.individualColonies
            if options.normalize
                warning('individual colonies not compatible with normalize')
            end
            for pi = 1:numel(positions)
                for cii = 1:numel(options.nucChannels)
                    ci = options.nucChannels(cii) + 1;
                    plot(r(1:end-binmarg), positions(pi).radialProfile.NucAvgSeg(1:end-binmarg,ci));%, 'Color', options.colors(cii,:)); 
                end
            end
        end
    end
        
    if ~isempty(options.cytChannels)
        
        for cii = 1:numel(options.cytChannels)
            ci = options.cytChannels(cii) + 1;
            if options.std
                errorbar(r, cyt_profile(:,ci), cyt_profile_std(:,ci),'--','LineWidth',2, 'Color', options.colors(cii,:)); 
            else
                plot(r, cyt_profile(:,ci),'--','LineWidth',2, 'Color', options.colors(cii,:))
            end
            legendentries = [legendentries, ['cyt ' meta.channelLabel{ci}]];
        end
    end
    hold off
    
    if options.legend
        legend(legendentries, 'Location','NorthEast');
    end

    % make it pretty
    fgc = 'k';
    bgc = 'w';
    graphbgc = 1*[1 1 1]; 

    xlim([0 positions(1).radiusMicron]);
    xlabel('edge distance ( um )', 'FontSize',fs,'FontWeight','Bold','Color',fgc)
    ylabel('intensity (a.u.)', 'FontSize',fs,'FontWeight','Bold','Color',fgc);

    set(gcf,'color',bgc);
    set(gca, 'LineWidth', 2);
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
    set(gca,'XColor',fgc);
    set(gca,'YColor',fgc);
    set(gca,'Color',graphbgc);
    
    res = struct('r',r, 'nuc_profile',nuc_profile, 'nuc_profile_std',nuc_profile_std,...
        'cyt_profile', cyt_profile, 'cyt_profile_std', cyt_profile_std);
end