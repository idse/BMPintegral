function [newCellData, chain] = combineZtest2(cellData, meta, opts)
% create new cellData merging info from different zslices

%options
%IoU = intersection over union threshold (sort of)
IoU = opts.IoU;
%convert microns to pixels or frames for other options
maxZslices = floor(opts.maxZsize/meta.zres);
maxZslices = max([maxZslices,2]);
fprintf('Max z slices = %d slices = %g um\n',maxZslices,maxZslices*meta.zres)
maxCentroidDistance = opts.maxCentroidDistance/meta.xres;
fprintf('maxCentroidDistance = %g um = %g pixels\n',...
    maxCentroidDistance*meta.xres,maxCentroidDistance)

% match nuclei between frames
nucChannel = meta.nucChannel;
nz = meta.nZslices;

% determine first non-empty slice; nei = non-empty index
nei = find(cellfun(@(x) ~isempty(x), {cellData.XY}),1,'first');
fprintf('First non-empty z index = %d\n',nei)

%number of cells in the first frame; use this to keep track of the number
%of uniquely identified nuclei across slices
cn = size(cellData(nei).XY,1);
chain = cell(cn,1);
for jj = 1:cn
    %chain{jj} = [z slice, cell index, nuclear marker level, nuclear area]
    nucLevel = cellData(nei).nucLevel(jj,meta.nucChannel+1);
    nucArea = cellData(nei).nucArea(jj);
    chain{jj} = [nei,jj,nucLevel,nucArea];
end
%for each slice, keep track of the chain to which each cell is assigned;
%for the first slice, each cell initializes its own chain
cidxs = cell(nz,1);
cidxs{nei} = (1:cn)';

%do the assignment as a linear assignment problem between adjacent z-slices
%analogously to frame-frame linking for the single cell tracking problem,
%but costs are based on how closely objects overlap
for ii = nei:nz-1
    XY1 = cellData(ii).XY;
    XY2 = cellData(ii+1).XY;
    %number of cells in each frame
    l1 = size(XY1,1); l2 = size(XY2,1); %(dont confuse l1 for 11 or l2 for 12)
    %don't try to link if one or both frames have no cells
    if (l1 > 0) && (l2 > 0)
        % impose minimal nubmer of cells (crash if ncells < k)
        k = min([l1,l2,3]);
        A11 = Inf(l1,l2);
        
        %only calculate overlap for a few possible candidates for each nucleus
        %but do it both ways (neighbors in slice ii+1 of a cell in slice ii and
        %neighbors in slice ii of a cell in slice ii+1)
        [Idx, D] = knnsearch(XY1,XY2,'K',k);
        for jj = 1:l2
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
        for jj = 1:l1
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
        A12 = Inf(l1,l1);
        A12(eye(l1)==1) = 1/IoU;
        A21 = Inf(l2,l2);
        A21(eye(l2)==1) = 1/IoU;

        %final cost matrix
        CM = [...
            A11, A12;...
            A21, A11'];
        %do the optimization
        [CM_indices, ~] = lapjv(CM);
        %parse resulting assignments
        target_indices = CM_indices(1:l1);
        source_indices = CM_indices((l1+1):(l1+l2)) - l2;
        %link rejection is denoted by index of -1
        target_indices(target_indices > l2) = -1;
        source_indices(source_indices < 1) = -1;

        %cells in frame ii+1 that were assigned to a cell in frame ii get added
        %to the chain that that cell is a part of
        C = NaN(l2,1);
        for jj = 1:l1
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
        cidxs{ii+1} = (cn+1):(cn+l2);
        for jj = 1:l2
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

% reduce cellData to have info for only one z slice for each cell
%make each nucleus have one set of pixels, defined in the slice in which
%its average intensity is highest
nCells = length(chain);
newCellData = struct;
fields = {'cytMaskArea','cytLevel','NCratio','XY','nucArea',...
    'nucOrientation','nucMajorAxis','nucMinorAxis','nucCircularity','nucZ'};
%initialize each field
for fi = 1:length(fields)
    if isfield(cellData(nei),fields{fi})
        newCellData.(fields{fi}) = zeros(nCells,size(cellData(nei).(fields{fi}),2));
    end
end
fields = fieldnames(newCellData);

%for each chain, find the slice for which the nucleus is brightest in the 
%nuclear channel and use data from that slice for the cell
for jj = 1:nCells
    cellinfo = chain{jj};
    [~, I] = max(cellinfo(:,3));
    zi = cellinfo(I,1);
    ci = cellinfo(I,2);

    for fi = 1:length(fields)
        newCellData.(fields{fi})(jj,:) = cellData(zi).(fields{fi})(ci,:);
    end
    
    %collect nuclear level in each channel and z slice for this cell in an
    %array
    nucLevel = zeros(size(cellinfo,1),size(cellData(zi).nucLevel,2));
    cytLevel = nucLevel;
    NCratio = nucLevel;
    for zi = 1:size(cellinfo,1)
        z = cellinfo(zi,1);
        ci = cellinfo(zi,2);
        nucLevel(zi,:) = cellData(z).nucLevel(ci,:);
        cytLevel(zi,:) = cellData(z).cytLevel(ci,:);
        NCratio(zi,:) = cellData(z).NCratio(ci,:);
    end
    %for each channel, keep the maximum intensity across z slices
    [~,I] = max(nucLevel,[],1);
    for ci = 1:size(nucLevel,2)
        newCellData.nucLevel(jj,ci) = nucLevel(I(ci),ci);
        newCellData.cytLevel(jj,ci) = cytLevel(I(ci),ci);
        newCellData.NCratio(jj,ci) = NCratio(I(ci),ci);
    end
end

background = cell2mat({cellData.background}');
newCellData.background = mean(background,1);



end