function dataCombined = combineQPCR(dataDir, filenames, tolerance)
    % combineQPCR(dataDir, filenames, tolerance)
    % tolerance: max std
    
    % SOP : StepOne Plus
    data = {};
    targetsComb = {};
    samplesComb = {};
    
    for i = 1:numel(filenames)
        
        fullfile(dataDir,filenames{i})
        data{i} = readSOPdata3(fullfile(dataDir,filenames{i}), tolerance);
        
        targetsComb = unique(cat(1,targetsComb,data{i}.targets),'stable');
        samplesComb = unique(cat(1,samplesComb,data{i}.samples),'stable');
    end

    %data{1}.targets{3} = 'NANOG';

    % collect CT values etc
    CTmean = [];
    CTstd = [];
    Ntargets = 0;
    Nsamples = 0;
    targets = [];
    samples = [];
    
    CTmeanRaw = [];
    CTstdRaw = [];
    outliers = [];
    highstd = [];
    
    if numel(targetsComb) > numel(data{1}.targets) &&...
        numel(samplesComb) > numel(data{1}.samples)

        error('either samples or targets have to match between plates');
        
    elseif numel(targetsComb) >= numel(data{1}.targets) &&...
            numel(samplesComb) == numel(data{1}.samples)
        
        if numel(filenames) > 1
            disp('combining multiple plates with different targets, same samples');
        end
        for i = 1:numel(data)

            D = data{i};

            CTmean = cat(2, CTmean, D.CTmean);
            CTstd = cat(2, CTstd, D.CTstd);

            Ntargets = Ntargets + D.Ntargets;
            targets = cat(1, targets, D.targets);

            CTmeanRaw = cat(2, CTmeanRaw, D.CTmeanRaw);
            CTstdRaw = cat(2, CTstdRaw, D.CTstdRaw);
            outliers = cat(2,outliers, D.outliers);
            highstd = cat(2,highstd, D.highstd);
            
            if i == 1
                Nsamples = D.Nsamples;
                samples = D.samples;
            elseif Nsamples ~= D.Nsamples
                error('samples dont match');
            end
        end

    elseif numel(targetsComb) == numel(data{1}.targets) &&...
            numel(samplesComb) > numel(data{1}.samples)

        disp('combining multiple plates with same targets, different samples');

        for i = 1:numel(data)

            D = data{i};

            CTmean = cat(1, CTmean, D.CTmean);
            CTstd = cat(1, CTstd, D.CTstd);

            Nsamples = Nsamples + D.Nsamples;
            samples = cat(1, samples, D.samples);

            CTmeanRaw = cat(1, CTmeanRaw, D.CTmeanRaw);
            CTstdRaw = cat(1, CTstdRaw, D.CTstdRaw);
            outliers = cat(1,outliers, D.outliers);
            highstd = cat(1,highstd, D.highstd);
            
            if i == 1
                Ntargets = D.Ntargets;
                targets = D.targets;
            elseif Ntargets ~= D.Ntargets
                error('samples dont match');
            end
        end
    end

    dataCombined = struct(...
        'targets',{targets},'samples',{samples},'CTmean',CTmean,'CTstd',CTstd,...
        'Nsamples', Nsamples, 'Ntargets', Ntargets,...
        'outliers', outliers, 'highstd', highstd,...
        'CTmeanRaw', CTmeanRaw, 'CTstdRaw', CTstdRaw);
    
    % to avoid misinterpretation of missing samples
    dataCombined.CTmean(dataCombined.CTmean == 0) = NaN;
end