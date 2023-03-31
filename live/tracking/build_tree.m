function trees = build_tree(tracks, splits)
%Given a collection of tracks and an array with information about which 
%tracks split from which others, build a tree with a collection of compound 
%tracks

%tracks is a cell array where each entry is a ntimepoints x 1 vector
%containing the index of the particle in that track at each time, with NaNs
%whenever the track is not present

%splits is an nx3 vector where each row gives the split track, the track it
%split from, and the time of the split - [split_track, mid_track, time]

%% 
%find tracks that are not split from other tracks:
n_tracks = size(tracks,1);
roots = transpose(1:n_tracks);
roots = roots(~ismember(roots, splits(:,1)));
%make a tree object that holds all the compound tracks; each compound track
%begins at a node connected to the root of the tree
trees = tree([]);
%Make a map variable that gives the correspondence between nodes of the
%tree and indices of tracks; the root node does not correspond to any
%track, so zero (arbitrary) is mapped to one (index of root)
maps = containers.Map(0,1);
%Add the content of all roots of compound tracks; the node index for each
%of these is one minus the index of the track in new_tracks
for k = 1:length(roots)
    trees = trees.addnode(1, tracks{roots(k)});
    %in the map for each point in roots, map the track index 
    maps(roots(k)) = nnodes(trees);
end

for i = 1:size(roots,1)
    criterion = true;
    %parent gives track index of parent track
    parent = roots(i);
    children = splits(splits(:,2) == parent, 1);
    %Make a list with a column of indices of tracks splitting off on the
    %left and indices of the tracks they split from on the right
    clist = [children, parent*ones(size(children,1),1)];
    loop_count = 0;
    while criterion
        plist = [];
        for j = 1:size(clist,1)
            idx = clist(j,1); %index of current track
            idxParent = clist(j,2); %index of parent track
            parentNode = maps(idxParent); %node corresponding to parent track

            %Add a child to parent node corresponding to track at idxParent
            trees = trees.addnode(parentNode, tracks{idx});

            %The new node is at the last index in the tree, so its index is
            %given by the total # of nodes immediately after it is added
            newnode = nnodes(trees);

            %Add an entry to the map showing that the track corresponding
            %to index in sim_tracks is in newnode in the tree
            maps(idx) = newnode;
            %trees{i,3}(newnode) = idx;

            %Find children of current track and add them to plist beside
            %current node (building new clist)
            children = splits(splits(:,2) == idx,1);
            plist = vertcat(plist,[children,idx*ones(size(children,1),1)]);
        end
        %Update clist
        clist = plist;
        loop_count = loop_count + 1;
        if loop_count > 100
            error('incorrect track assignment')
        end

        %If none of the tracks from this layer had children, then we are
        %done with this tree
        if isempty(clist)
            criterion = false;
        end
    end
end


end