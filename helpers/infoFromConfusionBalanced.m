function [MI, accuracy] = infoFromConfusionBalanced(M)
%estimate mutual information given a confusion matrix of true and assigned
%labels in a classification problem
%assume M(i,j) gives the number of samples of true class i with assigned
%class j

%balance the classes by normalizing each row to itself
m = size(M,1);
for ii = 1:m
    M(ii,:) = M(ii,:)/sum(M(ii,:));
end
%normalize the entire matrix so that entries are probabilities
M = M/sum(M,'all');

%accuracy is the trace of the normalized matrix
accuracy = trace(M);

%mutual information is found as sum(P(u,uhat)*log2(P(u,uhat)/(P(u)*P(uhat)))
MI = 0;
for ii = 1:m
    for jj = 1:m
        %increment
        inc = M(ii,jj)*log2(M(ii,jj)/(sum(M(:,jj))*sum(M(ii,:))));
        %if M(ii,jj) = 0, inc will be NaN because logarithms of 0 are
        %undefined, but we use the convention that 0*log(0) = 0
        if ~isnan(inc)
            MI = MI + inc;
        end
    end
end

end