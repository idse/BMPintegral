function M = confusionFromLabels(y_true,y_test)

%category labels
cat = unique([y_true(:); y_test(:)]);
%number of input categories
q = length(cat);
%initialize confusion matrix
M = zeros(q,q);
%fill the confusion matrix M; %M(i,j) gives the fraction of inputs from
%category i that are decoded as category j
for ii = 1:q
    for jj = 1:q
        M(ii,jj) = sum(y_true == cat(ii) & y_test == cat(jj));
    end
end


end