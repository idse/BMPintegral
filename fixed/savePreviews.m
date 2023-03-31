function savePreviews(images, experimentMeta, IlimAll, postfix)

P = experimentMeta;

if ~exist('postfix','var')
    postfix = [];
end

Nrows = size(images,1); %P.conditions
Ncols = size(images,2);
screensize = get( 0, 'Screensize' );

for overview = [true false]

if overview
    figure('Position', screensize);
end

for i = 1:Nrows
    for j = 1:Ncols
        if ~isempty(images{i,j})
            for k = 1:numel(images{i,j})

                li = P.channelLabelIdx{i,j};
                
                MIP = max(images{i,j}{k},[],4);
                for ci = 1:size(MIP,3)
                   MIP(:,:,ci) = imadjust(MIP(:,:,ci),IlimAll{li}(:,ci));
                end

                MIPRGB = zeros([size(MIP,1) size(MIP,2) 3],'uint16');
                for ci = P.channels
                    if P.channelColor{li}(ci) > 0 && P.channelColor{li}(ci) < 4
                        MIPRGB(:,:,P.channelColor{li}(ci)) = MIPRGB(:,:,P.channelColor{li}(ci)) + MIP(:,:,ci);
                    elseif P.channelColor{li}(ci) == 0
                        for cii = 1:3
                            MIPRGB(:,:,cii) = MIPRGB(:,:,cii) + MIP(:,:,ci);
                        end
                    end
                end

                ysize = size(MIP,1);
                %---
                channelLabels = P.channelLabels{li};
                colors = [[1 1 1]*0.9; 1 0.5 0.5; 0.5 1 0.5; 0.5 0.5 1];
                
                if overview
                    subplot_tight(Nrows,Ncols,sub2ind([Ncols Nrows],j,i));
                else
                    figure,
                end
                
                if overview
                    fs = 5;
                else
                    fs = 20;
                end
                
                imshow(MIPRGB)
                for ci = P.channels
                    text(50 + (ci-1)*200, ysize-75, channelLabels{ci},...
                        'FontSize',fs,'Color',colors(P.channelColor{li}(ci)+1,:),'FontWeight','bold')
                end
                text(50,50, P.conditions{i,j}{k},'FontSize',fs,'Color',[1 1 1]*0.99,'FontWeight','bold','Interpreter','none');
                
                if ~overview
                    [~,barefname,~] = fileparts(P.fileNames{i,j}{k});
                    filename = [barefname '_preview' postfix '.jpg'];
                    export_fig(fullfile(P.dataDir,filename));%,'-m0.5' for smaller );%-native for bigger
                    close; 
                end
            end
        else
            if overview
                disp([num2str(Ncols) ' ' num2str(Nrows) 'empty']);
                subplot_tight(Nrows,Ncols,sub2ind([Ncols Nrows],j,i));
                imshow(ones(size(MIP,[1,2])));
                condition = P.conditions{i,j};
                if isempty(condition), condition = 'NA'; end
                text(50,50,condition,'FontSize',fs,'Color','k','FontWeight','bold','Interpreter','none');
            end
        end
    end
end
if overview
    export_fig(fullfile(P.dataDir,['overview' postfix '.jpg']));
    close;
end
end