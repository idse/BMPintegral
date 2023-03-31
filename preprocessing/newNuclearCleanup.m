function newseg = newNuclearCleanup(seg, fgChannels, divChannels, opts)

maxArea = opts.decompopts.maxArea;
tau1 = opts.decompopts.tau1;
tau2 = opts.decompopts.tau2;

segChannel = fgChannels(~ismember(fgChannels, divChannels));
nucseg = seg == segChannel;
divmask = ismember(seg, divChannels);

overlap = zeros(size(seg));
for ii = 1:length(fgChannels)
    overlap = overlap + bwmorph(seg == fgChannels(ii),'thicken',2);
end
overlap = overlap > 1;

%transfer large objects from divmask to nucseg
CC = bwconncomp(divmask);
sizes = cellfun(@numel, CC.PixelIdxList);
idxs = CC.PixelIdxList(sizes > maxArea);
idxs = cell2mat(idxs(:));
divmask(idxs) = 0;
nucseg(idxs) = 1;

newseg = separate_fused(nucseg, tau1, tau2, opts.decompopts);
newseg = (newseg | divmask) & ~overlap;
newseg = nuclearCleanup(newseg, opts.cleanupOptions);

%remove big objects
if isfield(opts.cleanupOptions,'maxArea')
    props = regionprops(newseg,'Area','PixelIdxList');
    idxs = cell2mat({props([props.Area]>opts.cleanupOptions.maxArea).PixelIdxList}');
    newseg(idxs) = 0;
end

end