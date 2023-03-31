function data = readSOPdata3(filename, tolerance)
    % tolerance: max std
    %
    % output:
    % struct data with fields
    % - targets:        cell array of target names
    % - samples:        cell array of sample names
    % - CTmean
    % - CTstd
    % - Nsamples
    % - Ntargets
    % - outliers        replicate that lies outside of tolerance
    % - highstd         if all replicates are bad, this is true
    % - CTmeanRaw       if one replicate is excluded this contains original
    %                   mean, while CTmean contain mean after exclusion
    % - CTstdRaw
    
    fid = fopen(filename,'r','n','UTF-8');

    tline = fgetl(fid);

    firstlinefound = 0;
    rawdata = [];
    q = 1;
    while(ischar(tline))

        if firstlinefound == 0
             if contains(tline,'[Results]')
                 firstlinefound = 1;
                 %rawdata = strsplit(tline,'\t');
                 %disp(tline);
             end
        else
            if ~isempty(tline)
                x = strsplit(tline,'\t','CollapseDelimiters',false);
                % following bit serves to deal with data that has different
                % number of columns in different rows (has happened)
                if ~isempty(rawdata)
                    L = min(size(rawdata,2), numel(x));
                    rawdata = cat(1, rawdata(:,1:L), x(1:L));
                else
                    rawdata = x;
                end
                q = q + 1;
            else
                break;
            end
        end
        tline = fgetl(fid);
        tline = strrep(tline,',',', ');
    end
   
    Scoli = find(~cellfun(@isempty,strfind(rawdata(1,:),'Sample Name')),3,'first');
    Tcoli = find(~cellfun(@isempty,strfind(rawdata(1,:),'Target Name')),3,'first');
    Ccoli = find(~cellfun(@isempty,strfind(rawdata(1,:),'C')),3,'first');

    samples = unique(rawdata(2:end,Scoli),'stable');
    targets = unique(rawdata(2:end,Tcoli),'stable');

    Nsamples = numel(samples);
    Ntargets = numel(targets);

    % collect CT values

    CTmean = zeros([Nsamples Ntargets]);
    CTstd = zeros([Nsamples Ntargets]);
    
    outliers = false([Nsamples Ntargets]);
    highstd = false([Nsamples Ntargets]);
    CTmeanRaw = zeros([Nsamples Ntargets]);
    CTstdRaw = zeros([Nsamples Ntargets]);
    
	for si = 1:Nsamples
        for ti = 1:Ntargets

            ind = strcmp(rawdata(:,Scoli), samples{si}) & strcmp(rawdata(:,Tcoli), targets{ti});
            if any(ind)
                %disp([samples{si} '; ' targets{ti}]);
                CT = str2double(rawdata(ind,Ccoli(1)));
                CTstd(si,ti) = std(CT);
                if CTstd(si,ti) < tolerance
                    CTmean(si,ti) = mean(CT);
                else
                    M = distmat(CT);
                    good = ~all(M == 0 | M > tolerance | isnan(M));
                    if ~any(good)
                        disp(['std too high for: ' samples{si} ', ' targets{ti}]);
                        disp(['CT: ' num2str(CT',3)]);
                        
                        CTmean(si,ti) = nanmean(CT);
                        CTstd(si,ti) = nanstd(CT);
                        highstd(si,ti) = true;
                        
                    else
                        disp(['excluding outlier for: ' samples{si} ', ' targets{ti}]);
                        disp(['CT: ' num2str(CT',3)]);
                        
                        CTmean(si,ti) = mean(CT(good));
                        CTstd(si,ti) = std(CT(good));
                        
                        outliers(si, ti) = true;
                        CTmeanRaw(si,ti) = nanmean(CT);
                        CTstdRaw(si,ti) = nanstd(CT);
                    end
                end
            end
        end
    end

    data = struct('rawdata',{rawdata},...
        'targets',{targets},'samples',{samples},'CTmean',CTmean,'CTstd',CTstd,...
        'Nsamples', Nsamples, 'Ntargets', Ntargets,...
        'outliers', outliers, 'highstd', highstd,...
        'CTmeanRaw', CTmeanRaw, 'CTstdRaw', CTstdRaw);
end
