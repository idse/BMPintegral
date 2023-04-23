function [newedges, badedges, nodes, costs] = newgraphmergelinker(G,...
    points, areas, Divs, max_gap_time)
%Function to resolve merging events in single-cell tracks, since these are
%just artifacts of nuclei overlapping in some frames but not others
%If a merge is followed soon after (within the gap time cutoff) by a split,
%and that split does not appear to be due to cell division, assign each of
%the merge's input cells to one of the split's output cells, and remove
%intermediary links
%if there is not a split soon after, change the merge to a disappearance,
%assigning only one of the input cells to the output (just remove one of
%the links, leave others in place and make no additions)

%Initialize variables:
Din = indegree(G); %number of input edges to each node
Dout = outdegree(G); %number of output edges from each node
mergepoints = find(Din == 2); %nodes with two input edges
%collect area, division information in tall vectors (one row per node in
%tracking graph, so we don't have to index back and forth with time,id)
allpoints = cell2mat(points);
allA = cell2mat(areas);
allDiv = cell2mat(Divs);

%initialize these guys
newedges = [];
badedges = [];
costs = [];
nodes = [];

%for each merge, advance through successor points until (1) the
%max_gap_time is reached, (2) another merge is reached, (3) a split is
%reached, or (4) there are no more links (track end). If a split is reached
%and it is a non-dividing cell, match the two inputs to the two outputs;
%otherwise, match one and discard the other link
for ii = 1:length(mergepoints)
    idx = mergepoints(ii);
    ins = predecessors(G, idx);
    xyin = allpoints(ins,:);
    Ain = allA(ins);
    outs = idx;
    oldedges = zeros(max_gap_time, 2);
    
    %keep building/looping for the current track while criterion = true
    criterion = true;
    count = 0;
    %change flag to true if a non-dividing split is found
    flag = false;
    if Dout(idx) == 2
        criterion = false;
        if ~allDiv(idx)
            outs = successors(G,idx);
            flag = true;
        end
    end
    while criterion
        newidx = successors(G,idx);
        count = count + 1;
        
        if isempty(newidx)
            criterion = false;
            count = count - 1;
        else
            oldedges(count,:) = [idx, newidx];
            idx = newidx;
            if count == max_gap_time || Din(idx) == 2
                criterion = false;
            elseif Dout(idx) == 2
                criterion = false;
                if ~allDiv(idx)
                    outs = successors(G,idx);
                    flag = true;
                end
            end
        end
    end
    oldedges = oldedges(1:count,:);
    if any(oldedges == 0,'all')
        error('why is there still 0')
    end
    
    if flag
        %if we are linking merge inputs to split outputs, get rid of all
        %intermediary links (links in the merge, links between cells before
        %the split, links to the split) and add two links - one from each
        %input to one output
        newbadedges = [ins(1) mergepoints(ii); ins(2) mergepoints(ii);...
            oldedges; idx outs(1); idx outs(2)];
        badedges = [badedges; newbadedges];
        xyout = allpoints(outs,:);
        Aout = allA(outs);
        %c = [c(in1,out1), c(in1,out2);
        %     c(in2,out1), c(in2,out2)]
        c = pdist2(xyin,xyout,'squaredeuclidean').*(1 +...
            2*abs(Ain - Aout')./(Ain + Aout'));
        if c(1,1) + c(2,2) < c(1,2) + c(2,1)
            newedges = [newedges; ins(1), outs(1); ins(2), outs(2)];
            costs = [costs; c(1,1); c(2,2)];
        else
            newedges = [newedges; ins(1), outs(2); ins(2), outs(1)];
            costs = [costs; c(1,2); c(2,1)];
        end
    else
        %if not linking inputs to splitting outputs, save nodeID of point
        %merged to so the merge can be removed
        nodes = [nodes; outs];
    end
end





end