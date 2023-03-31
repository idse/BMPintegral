function IlimAll = getIntensityLimits(images, experimentMeta, tol)

    P = experimentMeta;

    IlimAll = cell(1,numel(P.channelLabels));
    for li = 1:numel(P.channelLabels)
        Nchannels = numel(P.channelLabels{li});
        IlimAll{li} = zeros([2 Nchannels]);
    end

    for i = 1:size(images,1)
        for j = 1:size(images,2)
            for k = 1:numel(images{i,j})
                MIP = max(images{i,j}{k},[],4);
                if numel(tol) == 1
                    Ilim = stretchlim(MIP,tol);
                else
                % specify separate limits for each channel
                    Ilim = zeros([2 numel(tol)]);
                    if size(tol,1) == 2 && (size(tol,2) > 1 || Nchannels == 1)
                        %specify upper and lower limits for each channel
                        for ci = 1:numel(tol)
                            Ilim(:,ci) = stretchlim(MIP(:,:,ci),tol(:,ci));
                        end
                    else
                        %specify a single limit percentage for each channel
                        for ci = 1:numel(tol)
                            Ilim(:,ci) = stretchlim(MIP(:,:,ci),tol(ci));
                        end
                    end
                end
                IlimAll{P.channelLabelIdx{i,j}} = max(IlimAll{P.channelLabelIdx{i,j}}, Ilim);
            end
        end
    end

end