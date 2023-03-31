function H = knn_entropy(varargin)
%use k-nearest-neighbors algorithm to estimate the shannon entropy in a
%collection of samples of a random variable

X = varargin{1};
if nargin == 1
    k = 1;
    dist = 'euclidean';
elseif nargin == 2
    k = varargin{2};
    dist = 'euclidean';
elseif nargin == 3
    k = varargin{2};
    dist = varargin{3};
end

% number of samples
n = size(X,2);
%dimension of samples
d = size(X,1);

%volume of a unit sphere of dimension d; may need to compute this
%iteratively to avoid bad numerics
cd = pi^(d/2)*2^d/gamma(d/2 + 1);

D = squareform(pdist(X',dist));

knnsum = 0;
for ii = 1:n
    trace = X(:,ii);
    idxs = [1:ii-1,ii+1:n];
    Y = sort(D(ii,idxs));
    if isinf(log2(2*Y(k)))
        error('this is inf')
    end
    knnsum = knnsum + log2(1*Y(k));
%     disp(knnsum/ii)
end

H = -psi(k) + psi(n) + log2(cd) + knnsum*d/n;
% error('why')


end