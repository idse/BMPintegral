function subtree = buildsubtree(tree_in,node)
%Given some starting node on an input tree, build a new tree that has the
%input node as its root and all its children and their children etc. as in
%the input tree
%This is very similar to build_tree but with slightly different inputs

%make a new tree with the content at tree_in(node) as its root content
subtree = tree(get(tree_in, node));
%make a map to track correspondences between the nodes of the parent tree
%and the nodes of the subtree
maps = containers.Map(node, 1);

%iteratively add children/dependencies
criterion = true;
clist = getchildren(tree_in, node);
clist = [clist; node*ones(size(clist))];
while criterion
    tlist = [];
    for j = 1:size(clist,2)
        idx = clist(1,j);
        parent = clist(2,j);
        content = get(tree_in, idx);
        new_parent = maps(parent);
        subtree = addnode(subtree, new_parent, content);
        maps(idx) = nnodes(subtree);
        children = getchildren(tree_in, idx);
        tlist = horzcat(tlist, [children; idx*ones(size(children))]);
    end
    clist = tlist;
    if isempty(clist)
        criterion = false;
    end
end

end