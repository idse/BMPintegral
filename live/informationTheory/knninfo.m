function H = knninfo(varargin)
%function to estimate mutual information between a data matrix of time
%traces and the labels between those traces
%syntax:
    %val = knninfo(X, labels)
    %val = knninfo(X, labels, k)
    %val = knninfo(X, labels, k, dist)
    %val = knninfo(X, labels, k, dist, fu)
%output:
    %val = estimate for mutual information between time traces in X and
    %labels in labels
%inputs:
    %X = d x N data matrix with N time traces each of dimension/length d
    %labels = N x 1 (or 1 x N?) vector of labels for the traces in X with m
    %distinct labels
    %k = number of nearest neighbors for k-nearest-neighbors algorithm
    %dist = distance metric; if not specified, use euclidean distance by
    %default
    %fu = distribution of input signals -> probability of observing each
    %signal value in labels; if not specified, use a uniform distribution
    %by default (m x 1)

%mutual information I(X,U) = H(X) - H(X|U) = H(U) - H(U|X)

%parse inputs
X = varargin{1};
labels = varargin{2}(:);
if nargin == 2
    k = 1;
    dist = 'euclidean';
elseif nargin == 3
    k = varargin{3};
    dist = 'euclidean';
elseif nargin == 4 || nargin == 5 || nargin == 6
    k = varargin{3};
    dist = varargin{4};
elseif nargin > 6
    error('wrong number of inputs')
end

[d, N] = size(X);
if nargin == 6
    cats = varargin{6};
else
    cats = unique(labels);
end
m = length(cats);

if nargin < 5
    qs = (1/m)*ones(m,1);
else
    qs = varargin{5};
end

D = squareform(pdist(X',dist));

%empirically estimate probabilities of each input given number of
%trajectories for each input and the total number of trajectories
ns = zeros(m,1);
for ii = 1:m
    ns(ii) = sum(labels == cats(ii));
end
% qs = ns/N;

%calculate d-dimensional unit sphere volume for knn estimate algorithm
%this part is numerically wonky because the numerator and denominator both
%become massive as d becomes large (similar to binomial coefficients)
%-> better way to calculate?
Vd = pi^(d/2)/gamma(d/2 + 1);

%for each input trajectory, estimate the conditional probability of finding
%that trajectory given each input/label
indices = (1:N)';
fRU = zeros(N,m);
for jj = 1:m
    midxs = labels == cats(jj);
    for ii = 1:N
        %for the row of the distance matrix corresponding to trace ii, get
        %the elements corresponding to category jj; if trace ii is in
        %category jj, remove it from the list
        idxs = midxs & indices ~= ii;
        Y = sort(D(ii, idxs));
        %estimate conditional probability
        if length(Y) >= k && Y(k) > 0
            fRU(ii,jj) = k/(ns(jj)*Vd*(Y(k)^d));
        else
            fRU(ii,jj) = 1;
        end
    end
end

%use the conditional probability matrix and the label probabilities to find
%the marginal probability of each output trace (w/out conditioning on input
%value)
fR = fRU*qs;
qn = cell(1,m);
for jj = 1:m
    qn{jj} = repmat(qs(jj)/ns(jj),1,ns(jj));
end
qn = cell2mat(qn);
HR = -qn*log2(fR);
% HR = -1/N*sum(log2(fR));

HRU = 0;
for jj = 1:m
    idxs = find(labels == cats(jj));
    for ii = 1:length(idxs)
        idx = idxs(ii);
        HRU = HRU - qs(jj)/ns(jj)*log2(fRU(idx,jj));
    end
end

H = HR - HRU;

end