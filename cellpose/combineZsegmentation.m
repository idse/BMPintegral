function [newCellData, chain] = combineZsegmentation(cellData, meta, opts)
% create new cellData merging info from different zslices

%options
if ~isfield(opts,'Zmethod')
    opts.Zmethod = 'useNuclearMarker';
end
if ~isfield(opts,'minZsize')
    opts.minZsize = 1;
end
if ~strcmp(opts.Zmethod,'useNuclearMarker')
    warning('z option: using maximum nuclear intensity separately in each channel\n change if not using dragonfly water objective')
end
%IoU = intersection over union threshold (sort of)
IoU = opts.IoU;
%convert microns to pixels or frames for other options
maxZslices = max([floor(opts.maxZsize/meta.zres),1]);
minZslices = floor(opts.minZsize/meta.zres);
fprintf('Max z slices = %d slices = %g um\n',maxZslices,maxZslices*meta.zres)
fprintf('Min z slices = %d slices = %g um\n',minZslices,minZslices*meta.zres)
maxCentroidDistance = opts.maxCentroidDistance/meta.xres;
fprintf('max centroid distance = %g um = %g pixels\n',...
    maxCentroidDistance*meta.xres,maxCentroidDistance)

% match nuclei between frames
nucChannel = meta.nucChannel;
nz = meta.nZslices;

% determine first non-empty slice; nei = non-empty index
nei = find(cellfun(@(x) ~isempty(x), {cellData.XY}),1,'first');
fprintf('First non-empty z index = %d\n',nei)

%number of cells in the first frame; use this to keep track of the number
%of uniquely identified nuclei across slices
cn = size(cellData(nei).XY,1); %cn = cell number
%build chains of linked nuclear masks
chain = cell(cn,1);
for jj = 1:cn
    %chain{jj} = [z slice, cell index, nuclear marker level, nuclear area]
    nucLevel = cellData(nei).nucLevel(jj,meta.nucChannel+1);
    nucArea = cellData(nei).nucArea(jj);
    chain{jj} = [nei,jj,nucLevel,nucArea];
end
%for each slice, keep track of the chain to which each cell is assigned;
%for the first z slice, each cell initializes its own chain
cidxs = cell(nz,1);
cidxs{nei} = (1:cn)';

%do the assignment as a linear assignment problem between adjacent z-slices
%with costs based on how closely objects overlap
for ii = nei:nz-1
    XY1 = cellData(ii).XY;
    XY2 = cellData(ii+1).XY;
    %number of cells in each frame
    nc1 = size(XY1,1); nc2 = size(XY2,1); %(dont confuse l1 for 11 or l2 for 12)
    %don't try to link if one or both frames have no cells
    if (nc1 > 0) && (nc2 > 0)
        % impose minimal nubmer of cells (crash if ncells < k)?
        k = min([nc1,nc2,3]);
        A11 = Inf(nc1,nc2);
        
        %only calculate overlap for a few possible candidates for each nucleus
        %but do it both ways (neighbors in slice ii+1 of a cell in slice ii and
        %neighbors in slice ii of a cell in slice ii+1)
        [Idx, D] = knnsearch(XY1,XY2,'K',k);
        for jj = 1:nc2
            list2 = cellData(ii+1).PixelIdxList{jj};
            for ki = 1:k
                ll = Idx(jj,ki);
                if D(jj,ki) < maxCentroidDistance
                    list1 = cellData(ii).PixelIdxList{ll};
                    A11(ll,jj) =...
                        min(numel(list1),numel(list2))/numel(intersect(list1,list2));
                end
            end
        end

        [Idx, D] = knnsearch(XY2,XY1,'K',k);
        for jj = 1:nc1
            list1 = cellData(ii).PixelIdxList{jj};
            for ki = 1:k
                ll = Idx(jj,ki);
                if D(jj,ki) < maxCentroidDistance
                    list2 = cellData(ii+1).PixelIdxList{ll};
                    A11(jj,ll) =...
                        min(numel(list1),numel(list2))/numel(intersect(list1,list2));
                end
            end
        end

        %make alternative cost matrices to reject links using the IoU threshold
        A12 = Inf(nc1,nc1);
        A12(eye(nc1)==1) = 1/IoU; %Inf square matrix with alternative costs on the diagonal
        A21 = Inf(nc2,nc2);
        A21(eye(nc2)==1) = 1/IoU;

        %final cost matrix
        CM = [...
            A11, A12;...
            A21, A11'];
        %do the optimization
        [CM_indices, ~] = lapjv(CM);
        %parse resulting assignments
        target_indices = CM_indices(1:nc1);
        source_indices = CM_indices((nc1+1):(nc1+nc2)) - nc2;
        %link rejection is denoted by index of -1
        target_indices(target_indices > nc2) = -1;
        source_indices(source_indices < 1) = -1;

        %cells in frame ii+1 that were assigned to a cell in frame ii get added
        %to the chain that that cell is a part of
        C = NaN(nc2,1);
        for jj = 1:nc1
            tidx = target_indices(jj);
            cidx = cidxs{ii}(jj);
            if tidx > 0
                nucLevel = cellData(ii+1).nucLevel(tidx,nucChannel+1);
                nucArea = cellData(ii+1).nucArea(tidx);
                chain{cidx} = [chain{cidx};ii+1,tidx,nucLevel,nucArea];
                C(tidx) = cidx;
            end
        end

        %nuclei in frame ii+1 with no link to one in frame ii get new chains
        %instead of being assigned to existing ones
        newinds = find(source_indices < 0);
        for jj = 1:length(newinds)
            cn = cn + 1;
            tidx = newinds(jj);
            nucLevel = cellData(ii+1).nucLevel(tidx,nucChannel+1);
            nucArea = cellData(ii+1).nucArea(tidx);
            chain{cn} = [ii+1,tidx,nucLevel,nucArea];
            C(tidx) = cn;
        end
        cidxs{ii+1} = C;
        
    else
        %if there are no cells in frame ii or ii + 1, every cell in frame
        %ii+1 gets a new chain
        cidxs{ii+1} = (cn+1):(cn+nc2);
        for jj = 1:nc2
            cn = cn + 1;
            nucLevel = cellData(ii+1).nucLevel(jj,meta.nucChannel+1);
            nucArea = cellData(ii+1).nucArea(jj);
            chain{cn} = [ii+1,jj,nucLevel,nucArea];
        end
    end
end

%deal with chains with more than maxZslices components
I = find(cellfun(@(x) size(x,1), chain) > maxZslices);
while ~isempty(I)
    for ii = 1:length(I)
        areas = chain{I(ii)}(:,4);
        areaDiffs = areas(2:end) - areas(1:end-1);
        secondDiffs = zeros(size(areas));
        secondDiffs(2:end-1) = areaDiffs(2:end) - areaDiffs(1:end-1);
        filler = min(secondDiffs)-1;
        secondDiffs(1) = filler; secondDiffs(end) = filler;
        [~,Idx] = max(secondDiffs);
        c1 = chain{I(ii)}(1:Idx,:);
        c2 = chain{I(ii)}(Idx+1:end,:);
        chain{I(ii)} = c1;
        cn = cn + 1;
        chain{cn} = c2;
    end
    I = find(cellfun(@(x) size(x,1), chain) > maxZslices);
end
%filter out chains with less than minZslices components
chain = chain(cellfun(@(x) size(x,1), chain) >= minZslices);

% reduce cellData to have info for only one z slice for each cell
%make each nucleus have one set of pixels, defined in the slice in which
%its average intensity is highest
nCells = length(chain);
newCellData = struct;
% fields = {'cytMaskArea','cytLevel','NCratio','XY','nucArea',...
%     'nucOrientation','nucMajorAxis','nucMinorAxis','nucCircularity','nucZ'};
fields = fieldnames(cellData);
fields = fields(cellfun(@(x) ~any(strcmp(x,{'PixelIdxList','background'})),fields));
%initialize each field
for fi = 1:length(fields)
    if isfield(cellData(nei),fields{fi})
        newCellData.(fields{fi}) = zeros(nCells,size(cellData(nei).(fields{fi}),2));
    end
end
% fields = fieldnames(newCellData);

%for each chain, find the slice for which the nucleus is brightest in the 
%nuclear channel and use data from that slice for the cell
for jj = 1:nCells
    cellinfo = chain{jj};
    %cellinfo(ii,:) = [z slice, cell index, nuclear marker level, nuclear area]
    [~, I] = max(cellinfo(:,3)); %row of cellinfo with highest nuclear level
    zi = cellinfo(I,1); %z slice of brightest nuclear marker
    ci = cellinfo(I,2); %cell index in that z slice

    for fi = 1:length(fields)
        newCellData.(fields{fi})(jj,:) = cellData(zi).(fields{fi})(ci,:);
    end
    
    newCellData.nucZ(jj) = zi;
    
    if ~strcmp(opts.Zmethod,'useNuclearMarker')
        %need to update this
        
        %collect nuclear level in each channel and z slice for this cell in an
        %array
        nucLevel = zeros(size(cellinfo,1),size(cellData(zi).nucLevel,2));
        cytLevel = nucLevel; NCratio = nucLevel;
        for zi = 1:size(cellinfo,1)
            z = cellinfo(zi,1); % z slice
            ci = cellinfo(zi,2); %cell index
            nucLevel(zi,:) = cellData(z).nucLevel(ci,:);
            cytLevel(zi,:) = cellData(z).cytLevel(ci,:);
            NCratio(zi,:) = cellData(z).NCratio(ci,:);
        end
        [nucmax, I] = max(nucLevel,[],1);
        if strcmp(opts.Zmethod,'old')
            %for each channel, keep the maximum nuclear intensity across z slices
            newCellData.nucLevel(jj,:) = nucmax;
        elseif strcmp(opts.Zmethod,'new')
            %for each channel, keep the nuclear intensity, cytoplasmic
            %intensity, and N:C ratio in the slice at which the nuclear
            %intensity is highest for that channel
            for chan = 1:size(nucLevel,2)
                newCellData.nucLevel(jj,chan) = nucLevel(I(chan),chan);
                newCellData.cytLevel(jj,chan) = cytLevel(I(chan),chan);
                newCellData.NCratio(jj,chan) = NCratio(I(chan),chan);
            end
        end
    end
end

background = cell2mat({cellData.background}');
newCellData.background = mean(background,1);



end