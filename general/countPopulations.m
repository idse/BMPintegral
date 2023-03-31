function counts = countPopulations(positions, meta, stats, dataDir, combo, conditionsidx)
    % counts.positivefraction_combo = nPositions x 2 x 2 x 2
    % give fraction positive for each position 
    % the 2nd to 4th dimensions are positive vs negative for each marker
    % 1 = positive, 2 = negative
    % similar for counts.positivefractionavg_combo nWells x 2 x 2 x 2
    % averaged over positions of a conditions/well 

    positive = {};

    Ncells = zeros([meta.nPositions 1]);
    NcellsTotal = zeros([meta.nWells 1]);
    NcellsStd = zeros([meta.nWells 1]);
    
    positivefraction = zeros([meta.nPositions meta.nChannels]);
    positivefractionavg = zeros([meta.nWells meta.nChannels]);
    positivefractionstd = zeros([meta.nWells meta.nChannels]);

    positivemean = zeros([meta.nPositions meta.nChannels]);
    positivemedian = zeros([meta.nPositions meta.nChannels]);
    positivemax = zeros([meta.nPositions meta.nChannels]);
    
    positivefraction_combo = zeros([meta.nPositions 2 2 2]);
    positivefractionavg_combo = zeros([meta.nWells 2 2 2]);
    positivefractionstd_combo = zeros([meta.nWells 2 2 2]);
    
    intensitymean_combo = zeros([meta.nWells meta.nChannels 2 2 2]);
    
    for condi = 1:meta.nWells

        condPos = meta.posPerCondition*(condi-1)+1:meta.posPerCondition*condi;

        for pi = condPos
            nucLevel = positions(pi).cellData.nucLevel;
            Ncells(pi) = size(nucLevel,1);
            background = positions(pi).cellData.background;    
            nucLevel = nucLevel - background;
            for ci = 1:meta.nChannels
                positive{pi, ci} = nucLevel(:,ci) > stats.thresholds(ci);
                positivemean(pi,ci) = mean(nucLevel(positive{pi, ci},ci));
                positivemedian(pi,ci) = median(nucLevel(positive{pi, ci},ci));
                m = max(nucLevel(positive{pi, ci},ci));
                if isempty(m)
                    m = 0;
                end
                positivemax(pi,ci) = m;
                positivefraction(pi,ci) = sum(positive{pi, ci})/numel(positive{pi, ci});
            end

            pppidx = positive{pi, 2} & positive{pi, 3} & positive{pi, 4};
            ppnidx = positive{pi, 2} & positive{pi, 3} & ~positive{pi, 4};
            pnpidx = positive{pi, 2} & ~positive{pi, 3} & positive{pi, 4};
            pnnidx = positive{pi, 2} & ~positive{pi, 3} & ~positive{pi, 4};
            nppidx = ~positive{pi, 2} & positive{pi, 3} & positive{pi, 4};
            npnidx = ~positive{pi, 2} & positive{pi, 3} & ~positive{pi, 4};
            nnpidx = ~positive{pi, 2} & ~positive{pi, 3} & positive{pi, 4};
            nnnidx = ~positive{pi, 2} & ~positive{pi, 3} & ~positive{pi, 4};
            
            positivefraction_combo(pi, 1, 1, 1) = sum(pppidx)/numel(positive{pi, combo(1)});
            positivefraction_combo(pi, 1, 1, 2) = sum(ppnidx)/numel(positive{pi, combo(1)});
            positivefraction_combo(pi, 1, 2, 1) = sum(pnpidx)/numel(positive{pi, combo(1)});
            positivefraction_combo(pi, 1, 2, 2) = sum(pnnidx)/numel(positive{pi, combo(1)});
            positivefraction_combo(pi, 2, 1, 1) = sum(nppidx)/numel(positive{pi, combo(1)});
            positivefraction_combo(pi, 2, 1, 2) = sum(npnidx)/numel(positive{pi, combo(1)});
            positivefraction_combo(pi, 2, 2, 1) = sum(nnpidx)/numel(positive{pi, combo(1)});
            positivefraction_combo(pi, 2, 2, 2) = sum(nnnidx)/numel(positive{pi, combo(1)});
            
            for ci = 1:meta.nChannels
                intensitymean_combo(pi, ci, 1, 1, 1) = mean(nucLevel(pppidx,ci));%./nucLevel(pppidx,1)); % DAPI normalize
                intensitymean_combo(pi, ci, 1, 1, 2) = mean(nucLevel(ppnidx,ci));%./nucLevel(ppnidx,1));
                intensitymean_combo(pi, ci, 1, 2, 1) = mean(nucLevel(pnpidx,ci));%./nucLevel(pnpidx,1));
                intensitymean_combo(pi, ci, 1, 2, 2) = mean(nucLevel(pnnidx,ci));%./nucLevel(pnnidx,1));
                intensitymean_combo(pi, ci, 2, 1, 1) = mean(nucLevel(nppidx,ci));%./nucLevel(nppidx,1));
                intensitymean_combo(pi, ci, 2, 1, 2) = mean(nucLevel(npnidx,ci));%./nucLevel(npnidx,1));
                intensitymean_combo(pi, ci, 2, 2, 1) = mean(nucLevel(nnpidx,ci));%./nucLevel(nnpidx,1));
                intensitymean_combo(pi, ci, 2, 2, 2) = mean(nucLevel(nnnidx,ci));%./nucLevel(nnnidx,1));

% log transform
%                 intensitymean_combo(pi, ci, 1, 1, 1) = mean(log(1+(exp(1)-1)*nucLevel(pppidx,ci)/stats.thresholds(ci)));%./nucLevel(pppidx,1)); % DAPI normalize
%                 intensitymean_combo(pi, ci, 1, 1, 2) = mean(log(1+(exp(1)-1)*nucLevel(ppnidx,ci)/stats.thresholds(ci)));%./nucLevel(ppnidx,1));
%                 intensitymean_combo(pi, ci, 1, 2, 1) = mean(log(1+(exp(1)-1)*nucLevel(pnpidx,ci)/stats.thresholds(ci)));%./nucLevel(pnpidx,1));
%                 intensitymean_combo(pi, ci, 1, 2, 2) = mean(log(1+(exp(1)-1)*nucLevel(pnnidx,ci)/stats.thresholds(ci)));%./nucLevel(pnnidx,1));
%                 intensitymean_combo(pi, ci, 2, 1, 1) = mean(log(1+(exp(1)-1)*nucLevel(nppidx,ci)/stats.thresholds(ci)));%./nucLevel(nppidx,1));
%                 intensitymean_combo(pi, ci, 2, 1, 2) = mean(log(1+(exp(1)-1)*nucLevel(npnidx,ci)/stats.thresholds(ci)));%./nucLevel(npnidx,1));
%                 intensitymean_combo(pi, ci, 2, 2, 1) = mean(log(1+(exp(1)-1)*nucLevel(nnpidx,ci)/stats.thresholds(ci)));%./nucLevel(nnpidx,1));
%                 intensitymean_combo(pi, ci, 2, 2, 2) = mean(log(1+(exp(1)-1)*nucLevel(nnnidx,ci)/stats.thresholds(ci)));%./nucLevel(nnnidx,1));
            end
        end
    end

    for condi = 1:meta.nWells

        condPos = meta.posPerCondition*(condi-1)+1:meta.posPerCondition*condi;

        NcellsTotal(condi) = sum(Ncells(condPos));
        NcellsStd(condi) = std(Ncells(condPos));
        
        positivefractionavg(condi,:) = mean(positivefraction(condPos,:),1);
        positivefractionstd(condi,:) = std(positivefraction(condPos,:),1);

        positivefractionavg_combo(condi,:,:,:) = mean(positivefraction_combo(condPos,:,:,:),1);
        positivefractionstd_combo(condi,:,:,:) = std(positivefraction_combo(condPos,:,:,:),1);
    end

    counts = struct(    'Ncells',Ncells,...
                        'positivemean',positivemean,...
                        'positivemedian',positivemedian,...
                        'positivemax',positivemax,...
                        'positivefraction',positivefraction,...
                        'positivefractionavg',positivefractionavg,...
                        'positivefractionstd',positivefractionstd,...
                        'positivefraction_combo',positivefraction_combo,...
                        'positivefractionavg_combo',positivefractionavg_combo,...
                        'positivefractionstd_combo',positivefractionstd_combo,...
                        'intensitymean_combo',intensitymean_combo);

    % VISUALIZE
    fgc = 'k';
    bgc = 'w';
    fs = 24;
    cases = 1:5;
    
    for casei = cases
    
        if casei == 1 % FRACTIONS OF EACH MARKER
            
            fnameprefix = 'positiveFractions_bar';
            legendstr = meta.channelLabel(2:4);
            vals = positivefractionavg(conditionsidx,2:4)'*100;
            errs = positivefractionstd(conditionsidx,2:4)'*100;
            ylabelstr = '+% of all cells';
            titlestr = [];
            
        elseif casei == 2 % NUMBERS OF EACH MARKER
            
            fnameprefix = 'positive_bar';
            legendstr = meta.channelLabel(2:4);
            vals = (diag(NcellsTotal(conditionsidx))*positivefractionavg(conditionsidx,2:4)/meta.posPerCondition)';
            errs = (diag(NcellsTotal(conditionsidx))*positivefractionstd(conditionsidx,2:4)/meta.posPerCondition)';
            ylabelstr = '+cells / colony';
            titlestr = [];
        
        elseif casei == 3 || casei == 4 % COMBO=SUBSET FRACTION || NUMBERS
            
            if numel(combo) == 3
                
                combosub = zeros([4 3]);
                combosub(:,combo(1)-1) = 1; % shift because we're not including DAPI here so 2->1
                combosub(:,combo(2)-1) = [1 1 2 2]; % 1 = +, 2 = -
                combosub(:,combo(3)-1) = [1 2 1 2];

                vals = [];
                errs = [];
                for i = 1:4 % +-, --, ++, -+

                    comboi = sub2ind([2 2 2],combosub(i,1),combosub(i,2),combosub(i,3));
                    dvals = positivefractionavg_combo(conditionsidx, comboi)';
                    derrs = positivefractionstd_combo(conditionsidx, comboi)';
                    
                    if casei == 4
                        dvals = NcellsTotal(conditionsidx)'.*dvals/meta.posPerCondition;
                        derrs = NcellsTotal(conditionsidx)'.*derrs/meta.posPerCondition;
                    elseif casei == 3
                        dvals = dvals*100;
                        derrs = derrs*100;
                    end
                    
                    vals = cat(1,vals, dvals);
                    errs = cat(1,errs, derrs);
                end

                legendstr = {[meta.channelLabel{combo(2)} '+' meta.channelLabel{combo(3)} '+'],...
                    [meta.channelLabel{combo(2)} '+' meta.channelLabel{combo(3)} '-'],...
                    [meta.channelLabel{combo(2)} '-' meta.channelLabel{combo(3)} '+'],...
                    [meta.channelLabel{combo(2)} '-' meta.channelLabel{combo(3)} '-']};
                
            elseif numel(combo)==2

                legendstr = {[meta.channelLabel{combo(2)} '+'],...
                    [meta.channelLabel{combo(2)} '-']};
                error('implement this');
            end
            
            if casei == 3
                fnameprefix = ['positiveFractionsCombo_' num2str(combo) '_bar'];
                ylabelstr = '% of all cells';
                titlestr = ['breakdown of ' meta.channelLabel{combo(1)} '+'];
            
            elseif casei == 4
                fnameprefix = ['positiveCombo_' num2str(combo) '_bar'];
                ylabelstr = 'cells / colony';
                titlestr = ['breakdown of ' meta.channelLabel{combo(1)} '+'];
            end
        
        elseif casei == 5 % CELL NUMBERS
            fnameprefix = 'cellnumbers_bar';
            vals = NcellsTotal(conditionsidx)'/meta.posPerCondition;
            errs = NcellsStd(conditionsidx)';
            ylabelstr = 'total cells / colony';
            legendstr = [];
        end

        %errorbar(vals', errs', 'LineWidth',2);
        errorbar_groups(vals, errs,...
                    'bar_names', meta.conditions(conditionsidx),...
                    'bar_width',0.75,'errorbar_width',0.5,...
                    'optional_errorbar_arguments',{'LineStyle','none','Marker','none','LineWidth',2});
        ylim([0 max(vals(:))+max(errs(:))]);
        legend(legendstr,'Location','best','FontSize',fs-10);
        if ~isempty(titlestr)
            title(titlestr);
        end
        
        ylabel(ylabelstr);
        set(gcf,'color',bgc);
        set(gca, 'LineWidth', 2);
        set(gca,'FontSize', fs)
        set(gca,'FontWeight', 'bold')
        set(gca,'XColor',fgc);
        set(gca,'YColor',fgc);
        set(gca,'Color',bgc);
        axis square
        saveas(gcf, fullfile(dataDir, [fnameprefix num2str(conditionsidx) '.png'])); 
        close;
    end
    
end