function tracks = graphtracks(G, ntime)

Din = indegree(G);
Dout = outdegree(G);
Frames = G.Nodes.frame;
[bins, binsizes] = conncomp(G,'Type','weak');
ntracks = length(binsizes);

tracks = cell(ntracks,1);

for jj = 1:ntracks
    if mod(jj,round(ntracks/20)) == 0
        fprintf('.')
    end
    list = find(bins == jj);
    start = list(Din(list) == 0);
    ends = list(Dout(list) == 0);
    numcols = length(ends);

    I = NaN(ntime,numcols);

    for ii = 1:numcols
        idxs = shortestpath(G,start,ends(ii));
        I(Frames(idxs),ii) = idxs;
    end
    
    tracks{jj} = I;
end
fprintf('\n')

end