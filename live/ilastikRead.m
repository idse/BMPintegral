function seg = ilastikRead(fname)

seg = squeeze(h5read(fname,'/exported_data'));
if contains(fname, 'Simple Segmentation')
    seg = seg == 1;
elseif contains(fname, 'Probabilities')
    %take first ilastik class and threshold if you have probabilities
    seg = squeeze(seg(1,:,:,:) > 0.5);
elseif contains(fname, 'Object Predictions')
    %no extra step needed for object classification
else
    error('Unknown output type')
end

seg = permute(seg,[2,1,3]);

end