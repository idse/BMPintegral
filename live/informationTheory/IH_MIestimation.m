clear all; close all;
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);

%% getting started: Monte-Carlo integration of 1D Gauss

f = @(X) exp(-X(:,1).^2/2);

L = 5;
dx = 0.01;
N = round(2*L/dx);
x = linspace(-L,L,N);

IRiemann = 0;
for i = 1:N
    IRiemann = IRiemann + f(x(i))*dx;
end

N = 10^6;
x = (rand([N 1]) - 0.5)*2*L; 
IMC = 2*L*mean(f(x));

IRiemann
Iexact = sqrt(2*pi)
IMC


%% integrate N-D Gaussian

f = @(X) exp(-sum(X.^2,2)/2);

D = 2;
L = 5;
N = 10^5;
Iexact = sqrt(2*pi)^D;

Ntrials = 20;
error = zeros([Ntrials 1]);

tic 
for t = 1:Ntrials

    x = (rand([N D]) - 0.5)*2*L; 
    IMC = (2*L)^D*mean(f(x));
    error(t) = sqrt((Iexact - IMC).^2)/Iexact;
end
toc

IMC
Iexact
%mean(error)

%% entropy of two correlated Gaussians

r = 0.5;
Sr = [1 r; r 1];
Nd = sqrt(det(2*pi*Sr)); 
Px = @(X) exp(-X'*inv(Sr)*X/2)/Nd; % probability density
h = @(X) -Px(X).*log2(Px(X)); % entropy density

%h([0 0]')

D = 2;
L = 5;
N = 10^7;
e = exp(1);

Hexact = log2(det(2*pi*e*Sr))/2;

x = (rand([N D]) - 0.5)*2*L; 
% HMC = 0;
% for i = 1:N
%     HMC = HMC + h(x(i,:)');
% end
% HMC = (2*L)^D*HMC/N;
%error = sqrt((Iexact - IMC).^2)/Iexact;

% vectorized
Px = @(X) exp(-sum(X'.*(inv(Sr)*X'),1)'/2)/Nd;
h = @(X) -Px(X).*log2(Px(X));
hmc = h(x);
hmc(isnan(hmc))=0;
HMC = (2*L)^D*mean(hmc);

HMC
Hexact
%error

%% mutual information of two correlated Gaussian

% I(X1 ; X2) = H(X) - H(X1 | X2)

% H(X1 | X2) = -\int dx1 dx2 P(X1,X2) log P(X1,X2)/P(X2)
% P(X2) = \int dx1 P(x1, x2)

% it seems like a hassle to calculate P(X2) here, maybe I shouldn't bother

Iexact = -log2(1-r^2)/2;

%% entropy of sum of two Gaussians with zero mean

Nd = @(s) sqrt(2*pi*s^2); 
Pxs = @(X,s) exp(-X'*X/(2*s^2))/Nd(s); % probability density

s1 = 2;
s2 = 1;
s = s1;

Px = @(X) (Pxs(X,s1) + Pxs(X,s2))/2;
h = @(X) -Px(X)*log2(Px(X)); % entropy density

D = 1;
L = 5*s;
N = 10^5;

x = (rand([N D]) - 0.5)*2*L; 
HMC = 0;
for i = 1:N
    HMC = HMC + h(x(i,:)');
end
HMC = (2*L)^D*HMC/N;

HMC

% only exact answer if the variances are equal because it reduces to a
% single Gaussian
Hexact = log2(2*pi*e*s^2)/2

%%
%---------------------------------------------------------------------------------------------
% NOW FOR REAL DATA
%---------------------------------------------------------------------------------------------

% raw data location, modify if not the same as the location of this script
% dataDir = scriptPath; 
dataDir = '/Users/idse/data_tmp/Seth/170726_BMPdoseResponse/'; 

load(fullfile(dataDir,'Data.mat'));

% The Data.mat file has 87x298 matrices A, N, R, X, Xnew, NS, NSnew, where each column is a time trace and 298x1 vector labels. 
% A, N, and R are the nuclear area, level of nuclear marker, and ratio of major to minor axis for each cell (used to determine if a cell was dividing). 
% X has the raw nuclear:cytoplasmic ratio, Xnew has the peak-removed N:C ratio, 
% NS has the nuclear:(average cytoplasmic) ratio to try to reduce noise due to the cytoplasmic mask, and 
% NSnew is the peak-removed matrix for NS
% 'labels' has a label of 1, 2, 3, or 4 for each trace, which correspond to doses of 0.1, 0.3, 1, and 30, respectively

%% visualize data

% plot traces and means
colors = lines(4);
clf 
hold on
for dose = 1:nDoses
    plot(NSnew(:,labels==dose),'Color',colors(dose,:));
end
for dose = 1:nDoses
    plot(mean(NSnew(:,labels==dose),2),'Color','w','LineWidth',4);
    plot(mean(NSnew(:,labels==dose),2),'Color',colors(dose,:),'LineWidth',2);
end
hold off

%%
% plot error bars
clf
hold on
for dose = 1:nDoses
    plot(mean(NSnew(:,labels==dose),2),'Color',colors(dose,:),'LineWidth',2);
end
for dose = 1:nDoses
    errorbar(mean(NSnew(:,labels==dose),2),std(NSnew(:,labels==dose),1,2),'Color',colors(dose,:),'LineWidth',3);
end
hold off
legend('0.1','0.3','1','30');


%%
% split up traces for each dose and calculate the correlation matrix and
% mean

nDoses = max(labels(:));

NSsplit = {};
S = {};    % covariance matrix for each dose S(d)
m = {};    % mean trace of each dose m(d)
ev = {};   % eigenvalues for each dose lambda(d)
invS = {};

for dose = 1:nDoses
    
    NSsplit{dose} = NSnew(:,labels==dose);
    S{dose} = cov(NSsplit{dose}');  % rows are time, columns are traces
    m{dose} = mean(NSsplit{dose},2);
    ev{dose} = eig(S{dose});
    
    % calculate inverse for later
    invS{dose} = inv(S{dose});
end

%%
% since the covariance matrix is close to singular, 
% we should truncate it to have a well defined inverse

% [V,D] = eig(A) such that A*V = V*D so inv(V)*A*V = D and this should be
% orthogonal so inv(V) = V^T 
% and to transform to a basis of eigenvectors, all vectors (traces) v -> V^Tv
Sall = cov(NSnew');
[V,D] = eig(Sall);

% sort in descending order
[d,ind] = sort(diag(D),'descend');
V = V(ind,ind);
Sall = Sall(ind,ind);

%%
%--------------------------------------------------------------------------
% for sanity check first truncate the full covariance matrix
% and compute the entropy of the traces in a single Gaussian approximation
%--------------------------------------------------------------------------

% number of principle components to truncate at
nPC = 3;

% truncated covariance matrix in PC basis
Sall_n = V'*Sall*V;
Sall_n = Sall_n(1:nPC,1:nPC);
tracesall_n = V'*NSnew;
tracesall_n = tracesall_n(1:nPC,:);
mall_n = mean(tracesall_n,2);
invSall_n = inv(Sall_n);
    
% single Gaussian for all traces:
% normalization of Pxall
Nall = sqrt(det(2*pi*Sall_n));
% Pxall = @(X) exp(- (X-mall_n)'*invSall_n*(X-mall_n)/2)/Nall; 
% vectorized
Pxall = @(X) exp(-sum((X-mall_n).*(invSall_n*(X-mall_n)),1)/2)/Nall; 

% size of integration domain
L = 5*max(Sall_n(:));

% check normalization
N = 10^8;
x = (rand([dim N]) - 0.5)*2*L + mall_n;
tic
imc = Pxall(x);
imc(isnan(imc))=0;
IMC = (2*L)^dim*mean(imc)
toc

% define and compute entropy
hall = @(X) -Pxall(X).*log2(Pxall(X));
tic
imc = hall(x);
imc(isnan(imc))=0;
HMC = (2*L)^dim*mean(imc)
toc
Hexact_all = log2(det(2*pi*e*Sall_n))/2

%%
%--------------------------------------------------------------------------
% now truncate each dose separately 
% and compute the entropy of the traces in a sum of Gaussians approximation
%--------------------------------------------------------------------------

traces_n = {};
m_n = {};
S_n = {};
ev_n = {};
invS_n = {};

maxNPC = 10; % maximal number of principle components

HMC_n = zeros([maxNPC 1]); % entropy of Gaussian mixture
Hexact_all_n = zeros([maxNPC 1]); % entropy of combined single Gaussian
HXD_n = zeros([maxNPC 1]); % conditional entropy of Gaussian mixture
MI_n = zeros([maxNPC 1]); % mutual information

%-----------------------------------------------
% number of principle components to truncate at
for nPC = 1:maxNPC
%-----------------------------------------------

    HXD = 0;
    for dose = 1:nDoses

        % transformed and truncated traces
        traces_n{dose} = V'*NSsplit{dose};
        traces_n{dose} = traces_n{dose}(1:nPC,:);
        % mean of transformed traces
        m_n{dose} = mean(traces_n{dose},2);

        % transformed covariance matrices for each dose
        S_n{dose} = V'*S{dose}*V;
        S_n{dose} = S_n{dose}(1:nPC,1:nPC);

        ev_n{dose} = eig(S_n{dose});
        invS_n{dose} = inv(S_n{dose});
        
        % conditional entropy
        % H(X | d) = log(det(2pi e S(d)))/2
        HXd = log2(det(2*pi*e*S_n{dose}))/2;
        % H(X | D) = \sum_d H(x | d) P(d)
        HXD = HXD + HXd/nDoses;
    end
    HXD_n(nPC) = HXD;
    
    %----------------------------------------
    % single combined Gaussian for comparison 
    %----------------------------------------
    
    % truncated covariance matrix in PC basis
    Sall_n = V'*Sall*V;
    Sall_n = Sall_n(1:nPC,1:nPC);
    tracesall_n = V'*NSnew;
    tracesall_n = tracesall_n(1:nPC,:);
    mall_n = mean(tracesall_n,2);
    Hexact_all_n(nPC) = log2(det(2*pi*e*Sall_n))/2;
    
    %----------------------------------------
    % Gaussian mixture approximation
    %----------------------------------------

    % I(X | D) = H(X) - H(X | D)
    % H(X) = - sum_x P(x) log P(x) 
    % P(x) = 1/4 sum_d P(x|d) 

    % Gaussian for each dose
    % normalization:
    Nd = @(dose) sqrt(det(2*pi*S_n{dose})); 
    % Pxd = @(X,dose) exp(- (X-m_n{dose})'*invS_n{dose}*(X-m_n{dose})/2)/Nd(dose); 
    % vectorized:
    Pxd = @(X,dose) exp(-sum((X-m_n{dose}).*(invS_n{dose}*(X-m_n{dose})),1)/2)/Nd(dose); 
    % sum of Gaussians:
    Px = @(X) (Pxd(X,1) + Pxd(X,2) + Pxd(X,3) + Pxd(X,4))/4;

    % preserve integration domain of combined Gaussian
    L = 5*max(Sall_n(:));
    N = 10^8;
    dim = nPC;
    x = (rand([dim N]) - 0.5)*2*L + mall_n;
%     % check normalization
%     tic
%     imc = Px(x);
%     imc(isnan(imc))=0;
%     IMC = (2*L)^dim*mean(imc)
%     toc

    % entropy
    tic
    h = @(X) -Px(X).*log2(Px(X)); 
    imc = h(x);
    imc(isnan(imc))=0;
    HMC_n(nPC) = (2*L)^dim*mean(imc);
    toc
    
    MI_n(nPC) = HMC_n(nPC) - HXD_n(nPC);
end

% %% trying to find syntax to sum over arbitrary number of doses
% X = m{d};
% f = @(d) exp(- (X-m{d})'*inv(S{d})*(X-m{d}))/Nd;
% arrayfun(@f,1:nDoses,'uniformOutput',false); 

%% SAVE

%save('MIestimation_results');

%% plot entropy and mutual information as a function of principal components
nPC = 10;
plot(1:nPC, MI_n, '-x')
hold on
plot(1:nPC, HXD_n, '-x')
plot(1:nPC, HMC_n, '-x')
plot(1:nPC, Hexact_all_n, '-x')
xlabel('# principal components');
ylabel('information (bits)')
hold off
legend({'MI = H(X)-H(X|D)','H(X|D) (exact)','Gaussian mixture H(X) (MC)','combined entropy (exact)'},'Location','NorthWest')

saveas(gcf,'GaussianMI.png');

%% why does the exact Gaussian entropy decrease with increasing dimension?

% Artefact of the differential entropy:
% when covariances get smaller than 1/2*pi*e, the log becomes negative
% in H = log2(det(2*pi*e*S))/2 for the Gaussian.
%
% This is closely related to the fact that -P log P goes negative
% if P > 1, which only happens for probability densities, not
% probabilities. 
% For a Gaussian the maximal value of P is 1/sqrt(det(2*pi*S_n{dose}))
% so for 1D, P > 1 if sqrt(2*pi) sigma < 1, so sigma < 1/sqrt(2*pi) = 0.39
% that doesn't make the integral P log P negative yet (thats the 1/2*pi*e)
% but it is already a breakdown of the definition.
% Since the overall covariances are 2.62, 0.36, 0.12, things already break
% down for the second principal component.
% 
% so what you really need is sum P dx log (P dx) . I am not sure 
% what the continuum limit is but we could discretize all this.

%%
% plot principal components -> almost look like Fourier modes
colors = lines(maxNPC);
clf 
hold on
for n = 1:5%maxNPC
    plot(V(:,n),'Color',colors(n,:),'LineWidth',1);
end
hold off

%%
% let's project onto the principal components 
% -> same is true : d(mean(t)) ~ std(t) 

nPC = 3;
colors = lines(nDoses);
clf 
hold on
for dose = 1:nDoses
    trace_transform = V'*NSnew(:,labels==dose);
    trace_transform(nPC+1:end,:) = 0;
    trace_projected = V*trace_transform;
    errorbar(mean(trace_projected,2),std(trace_projected,1,2),'Color',colors(dose,:));
    [mean(mean(trace_projected,2)) mean(std(trace_projected,1,2))]
end
hold off

%% project on first PC, then calculate mean and std
% reasonable numbers come out
% so what is the problem?

PCproj = V(:,1)'*(NSnew - mean(NSnew,2));
for i = 1:4
    [mean(PCproj(:,labels==i)) std(PCproj(:,labels==i))]
end

% principle variance is still the same 
%S = cov((NSnew - mean(NSnew,2))');
S = cov(NSnew');
[V,D] = eig(S);

% but that's not the variance along the first principal component?
var(PCproj)
