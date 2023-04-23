function histories = graphSignalingHistories(P, G, fields, idxs)
%P is a Position object 
%G is a directed acyclic graph holding tracking links for live cell imaging
%data
%idxs gives the indices of cells in the final frame of the live data that
%were matched to cells in the in the fixed data; the order should be
%according to the order of the fixed cells (i.e. idxs(ii) = index of live
%cell matched to iith fixed cell); this can also be used to remove any
%other cells for which we don't want signaling histories, like cells that
%are dividing in the last frame, or to only make histories for those cells
%for which tracking has been manually validated

n = length(idxs);
nfields = length(fields);

%get the IDs of nodes in the graph from indices of cells in the last frame
nodeIDs = 1:G.numnodes;
nodeIDs = nodeIDs(end-(P.ncells(end)-1):end);
nodeIDs = nodeIDs(idxs);

cellData = struct;
for jj = 1:nfields
    %concatenate cell data for each field into a single matrix across time
    %points
    cellData.(fields{jj}) = cell2mat({P.cellData.(fields{jj})}');
end

%intialize a struct to hold signaling histories
histories(n) = struct;

for ii = 1:n
    %iterate backwards through the graph to build the signaling history
    %this is slow -> is there a better way using built-in functions?
    nodeID = nodeIDs(ii);
    cidxs = nodeID;
    criterion = true;
    while criterion
        prevID = predecessors(G, nodeID);
        
        if isempty(prevID)
            criterion = false;
        else
            nodeID = prevID;
            cidxs = [cidxs; nodeID];
        end
    end
    %reverse order
    cidxs = cidxs(end:-1:1);
    %add info for these nodes to the signaling history
    histories(ii).CellIdxs = G.Nodes.cellidx(cidxs);
    histories(ii).Time = G.Nodes.frame(cidxs);
    for jj = 1:nfields
        histories(ii).(fields{jj}) = cellData.(fields{jj})(cidxs,:);
    end
end




end