function [MI, accuracy] = infoFromConfusion(M)
%estimate mutual information given a confusion matrix of true and assigned
%labels in a classification problem

%normalize the matrix so that entries are probabilities; if they are
%already probabilities, sum(M,'all') = 1, so this will result in no change
M = M/sum(M,'all');
m = size(M,1);
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