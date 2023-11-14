clear; close all; clc

%% setup, load data
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
% dataDir = fullfile(scriptPath,'data');
dataDir = scriptPath;
%save figures in a subfolder of where the script is located
savedir = fullfile(scriptPath,'figures');
if ~exist(savedir,'dir'), mkdir(savedir); end

%tvec - time points at which SMAD4 and SOX2 are sampled
%SMAD4 - input SMAD4 dynamics
%SOX2 - resulting SOX2 dynamics
%ISL1 - end-point ISL1

load(fullfile(dataDir,'data','meta.mat'))
load(fullfile(dataDir,'data','SMAD4_SOX2_220905.mat'))
SOX2_220905 = SOX2;
SMAD4_220905 = SMAD4;
ntime = length(tvec);

load(fullfile(dataDir,'data','S4levels.mat'))
ncond = size(SOX2,2);

ISL1_220905 = load(fullfile(dataDir,'data','ISL1_220905.mat'));
ISL1_220830 = load(fullfile(dataDir,'data','ISL1_220830_combined.mat'));

SOX2_220830 = load(fullfile(dataDir,'data','SOX2_220830.mat'),'SOX2');
SOX2_220830 = SOX2_220830.SOX2;

SOX2avg = (SOX2_220830 + SOX2_220905)/2;
ISL1avg = (ISL1_220830.ISL1m + ISL1_220905.ISL1m)/2;
ISL1avgerr = (ISL1_220830.ISL1err + ISL1_220905.ISL1err)/2;

SMAD4_220629 = load(fullfile(dataDir,'data','SMAD4_220629'));
SMAD4_220629 = preprocess_SMAD4_220629(SMAD4_220629);
SMAD4_220629 = SMAD4_220629(1:ntime,:);
SMAD4_ideal = makeIdealS4(tvec,s4levels);

treatmentTime = 4;

fs = 28; lfs = 20; lw = 3;

%dset = 'avg', '220830', '220905'
dset = 'avg';
if strcmp(dset,'avg')
    SOX2 = SOX2avg;
    ISL1 = ISL1avg;
    ISL1err = ISL1avgerr;
elseif strcmp(dset,'220830')
    SOX2 = SOX2_220830;
    ISL1 = ISL1_220830.ISL1m;
    ISL1err = ISL1_220830.ISL1err;
elseif strcmp(dset,'220905')
    SOX2 = SOX2_220905;
    ISL1 = ISL1_220905.ISL1m;
    ISL1err = ISL1_220905.ISL1err;
end

% s4 = '220629';
s4 = 'idealS4';
if strcmp(s4,'idealS4')
    SMAD4 = SMAD4_ideal;
elseif strcmp(s4,'220629')
    SMAD4 = SMAD4_220629;
end
integrals = mean(SMAD4,1);

%set ISL1 to zero at mTeSR control
ISL1n = ISL1 - ISL1(end);
ISL1n(ISL1n < 0) = 0;

%initial SOX2 (should be 1 because of normalization)
y0 = mean(SOX2(tvec <= 0,:),'all');

%% early SOX2 slope vs SMAD4 integral
figpos = figurePosition(560,560);

tend = 16;
tmask = tvec > 0 & tvec < tend;
% levels = mean(SMAD4(tmask,:),1);
levels = mean(SMAD4(tmask,:),1)*tend/tvec(end);
xx = tvec(tmask); xx = xx(:);
xinv = pinv([xx,ones(length(xx),1)]);

figure
mbs = NaN(2,ncond);
for ii = 1:ncond
    yy = SOX2(tmask,ii);
    mb = xinv*yy;
    mbs(:,ii) = mb;
    cla
    plot(tvec,SOX2(:,ii),'o')
    hold on
    plot(tvec,mb(1)*tvec + mb(2))
    ylim([0,max(SOX2,[],'all')])%ylim([75,3500])
    xlim(tvec([1,end]))
    title(meta.conditions{ii})
    disp(meta.conditions{ii})
    fprintf('m = %g, b = %g\n',mb(1),mb(2))
    xline(16,'Color','k'); xline(0,'Color','k');
    % drawnow
    pause(0.25)
end

conds = {[1:8,17],9:17};
legstr = meta.conditions(conds{1});
legstr = strrep(legstr,'LDN',''); legstr = strrep(legstr,',42hr','');
figure('Position',figpos); hold on
for ii = 1:length(conds)
    plot(mbs(1,conds{ii}),'o','LineWidth',lw,'MarkerSize',15)
end
legend('42hr', '32hr','Location','northwest')
cleanSubplot(fs); axis square
xticks(1:length(conds{1})); xticklabels(legstr); xtickangle(45)
ylabel('SOX2 slope (au / hr)')
xlabel('BMPRi (nM)')

figure('Position',figpos); hold on
conds = {[1:8,17],9:17};
for ii = 1:length(conds)
    plot(levels(conds{ii}),tend*mbs(1,conds{ii}),'x','LineWidth',lw,'MarkerSize',15)
end

ll = levels';
linv = pinv([ll,ones(size(ll))]);
mb = tend*linv*mbs(1,:)';
ls = sort(ll);
plot(ls,mb(1)*ls + mb(2),'LineWidth',lw,'LineStyle','--','Color',[0,0,0,0.65])

legend('42hr', '32hr','Location','northeast')
cleanSubplot(fs); axis square

ylabel('\Delta SOX2 (au)') 
xlabel('SMAD4 (N:C) integral') 

savefigure(fullfile(dataDir,'figures',[dset,'_SOX2changevS4int']))

%% SMAD4-SOX2 inhibition function
%on short intervals along the SOX2 histories, find the average SMAD4
%signaling level, average SOX2 level, and fit a line to the SOX2 data on
%the interval to estimate the slope
%this allows us to directly estimate the form of the inhibition function of
%SMAD4 on SOX2
savedir = fullfile(dataDir,[dset,'_SOX2_',s4,'_SMAD4'],'SOX2');
if ~exist(savedir,'dir'), mkdir(savedir); end

margin = 0.1;

% conds = 1:8;
conds = 1:16;
nc = length(conds);

% pts = 2:40;
pts = 2:16;
window = 3;
npts = length(pts);
sox2avg = zeros(npts,nc);
smad4avg = zeros(npts,nc);
slopes = zeros(npts,nc);

for ci = 1:nc
    cidx = conds(ci);
    for ii = 1:npts
%         mask = tvec >= pts(ii) & tvec <= pts(ii+1);
        mask = tvec >= pts(ii) - window/2 & tvec <= pts(ii) + window/2;
        t = tvec(mask); sox2 = SOX2(mask,cidx); smad4 = SMAD4(mask,cidx);
        sox2avg(ii,cidx) = mean(sox2); smad4avg(ii,cidx) = mean(smad4);
        mb = [t',ones(size(t'))]\sox2;
        slopes(ii,cidx) = mb(1);
    end
end
x = smad4avg(:); y = sox2avg(:); z = slopes(:);

figure('Position',figurePosition(3*560,560))
subplot_tight(1,3,1,margin); hold on
colorscatter(x,y,z)
cleanSubplot; axis square; colormap turbo
xlabel('smad4'); ylabel('sox2'); h = colorbar; h.Label.String = 'dydt';
title('data')

%linear model
lfo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0,0],...
               'Upper',[Inf,Inf,Inf],...
               'StartPoint',[0.058, 0.058, 0.058]);
lft = fittype(@(beta,lambda,alpha,x,y) beta - lambda*x - alpha*y,...
    'independent', {'x', 'y'}, 'dependent', 'z', 'options', lfo);
[fitobject,gof] = fit([x,y],z,lft);

lambda = fitobject.lambda;
alpha = fitobject.alpha;
beta = fitobject.beta;

subplot_tight(1,3,2,margin); hold on
scatter3(x,y,z,'filled')

xx = linspace(0,1.3,500);
yy = linspace(0,1,500);
[XX,YY] = meshgrid(xx,yy);
surf(XX,YY,beta - lambda*XX - alpha*YY,'LineStyle','none')
xlabel('SMAD4'); ylabel('SOX2'); h = colorbar; h.Label.String = 'dydt';
cleanSubplot; axis square; colormap turbo
title('fit')
colorbar

subplot_tight(1,3,3,margin); hold on
zp = beta - lambda*x - alpha*y;
colorscatter(x,y,z - zp)
cleanSubplot; axis square; colormap turbo
xlabel('SMAD4'); ylabel('SOX2'); h = colorbar; h.Label.String = 'residuals';
title('residuals')

savefigure(fullfile(savedir,'planeFitTodSOX2dt.png'))

%% model for ISL1 accumulation in response to SMAD4 and SOX2
tmax = 42;
tspan = [0,tmax];
ic = 0;
%initialize parameters
alphai = 0.08;
lambdai = 0.12;

% opts = struct('betai',betai,'alphai',alphai,'Ksi',0.5,'Ks4',0.5,'n',2);
opts = struct(...
    'lambdai',  lambdai,...
    'alphai',   alphai,...
    'Ksi',      0.3,...
    'n',        4);
%fields over which to optimize
fields = {'lambdai','alphai','Ksi','n'};

nt = length(fields);

trainconds = 1:17;
nconds = length(trainconds);

%% fit ISL1 model parameters with real SOX2, SMAD4 inputs
%create a folder to save results, labeled with the current date
% c = clock;
% savedir = sprintf('%d%.2d%.2d',c(1),c(2),c(3));
% savedir = fullfile(dataDir,savedir(3:end),[s4,'_SMAD4'],'ISL1');
% savedir = fullfile(dataDir,[dset,'_SOX2_',s4,'_SMAD4'],'ISL1');
savedir = fullfile(dataDir,'ISL1_fit_n3');
if ~exist(savedir,'dir'), mkdir(savedir); end

%set hyperparameters
mu = zeros(1,nt); %zero vector for mean of gaussian rv
sigma = 5e-5; %variance determines the step size for each parameter 
% Sigma = diag(sigma*[1 1 1]); %diagonal covariance matrix
Sigma = diag(sigma*ones(1,nt)); %diagonal covariance matrix
tempmax = 0.1; %'temperature' determines the acceptance criterion
niters = 1000; %number of iterations per optimization run
nrepeats = 5; %number of independent repeats of the optimzation procedure

for ri = 1:nrepeats
    %randomly initialize parameters
%     theta = rand(nt,1);
    theta = opts2theta(opts,fields);
    newopts = theta2opts(theta,opts,fields);
    V = NaN(nconds,1);
    %simulate once to determine initial error
    parfor idx = 1:nconds
        cidx = trainconds(idx);
        S4 = SMAD4(:,cidx); S2 = SOX2(:,cidx); %#ok<PFBNS>
        [~,y] = ode45(@(t,y) ISL1_hybrid(t,y,tvec,S4,S2,newopts), tspan, ic); %#ok<PFTUSW>
        V(idx) = y(end);
    end
    
    errold = mean(sum((ISL1n(trainconds) - V).^2));
    fprintf('MSE = %g\n',errold)
    
    %initialize variables to store intermediate cost and parameter values
    thetas = zeros(nt,niters);
    costs = zeros(niters,1);
    oldcosts = zeros(niters,1);
    naccept = 0;
    close all
    tic
    %show the error over time
    % figure('Position',[1000 602 1250 750])
    figure('Position',figurePosition([1250 750]))
    for ii = 1:niters
        %update parameters with random step
        thetap = theta + mvnrnd(mu,Sigma,1)';
        %if any parameters stepped below zero, set them between zero and
        %their initial value instead
        halfp = 0.5*theta;
        thetap(thetap <= 0) = halfp(thetap <= 0);
        newopts = theta2opts(thetap,opts,fields);
        %run the simulation over treatment conditions in parallel
        V = NaN(nconds,1);
        parfor idx = 1:nconds
            cidx = trainconds(idx);
            S4 = SMAD4(:,cidx); S2 = SOX2(:,cidx); %#ok<PFBNS>
            [~,y] = ode45(@(t,y) ISL1_hybrid(t,y,tvec,S4,S2,newopts), tspan, ic); %#ok<PFTUSW>
            V(idx) = y(end);
        end
        %error with new parameter values
        errnew = mean(sum((ISL1n(trainconds) - V).^2),'all');
        
        costs(ii) = errnew;
        oldcosts(ii) = errold;

        %acceptance criterion
        temp = tempmax*(1 - ii/niters);
        if exp((errold - errnew)/temp) > rand(1)
            theta = thetap;
            errold = errnew;
            naccept = naccept + 1;
        end
        thetas(:,ii) = theta';
        %update the graph of costs every 10 iterations
        if mod(ii,10) == 0
            fprintf('.') %progress indicator
            clf
            semilogy(costs(1:ii),'LineWidth',2)
            hold on
            semilogy(oldcosts(1:ii),'LineWidth',2)
            hold off
            cleanSubplot
            xlabel('iterations'); ylabel('mean squared error')
            axis square
            drawnow
        end
        if mod(ii,500) == 0
            fprintf('\n')
        end
    end
    [cmin,I] = min(costs);
    theta = thetas(:,I);
    newopts = theta2opts(theta,opts,fields);

    fprintf('\nEfficiency = %g\n',naccept/niters)
    fprintf('MSE = %g\n',cmin)
    disp(newopts)
    toc
    
    %for each run, save a graph of proposed and accepted costs over
    %iterations and the associated parameter and cost values as mat files
    close all
    figure('Position',figurePosition([1250 750]))
    semilogy(costs,'LineWidth',2); hold on
    semilogy(oldcosts,'LineWidth',2); hold off
    xlabel('iterations'); ylabel('mean squared error')
    cleanSubplot; axis square
%     ylim([0,4])
    legend('proposed cost', 'accepted cost')
    saveas(gcf,fullfile(savedir,sprintf('costs_run%.2d.png',ri)))
    
    
    newopts = theta2opts(theta,opts,fields);
    save(fullfile(savedir,sprintf('paramValues_run%.2d',ri)),'newopts','thetas','fields','costs','oldcosts')
    
    margin = 0.02;
    close all
    % figure('Position',[10 550 2500 700])
    figure('Position',figurePosition([2500 700]))
    idx = 1;
    for ii = 1:nconds
        
        cidx = trainconds(ii);
        S4 = SMAD4(:,cidx); S2 = SOX2(:,cidx);
        [t,y] = ode45(@(t,y) ISL1_hybrid(t,y,tvec,S4,S2,newopts), tspan, ic);
%         [t,y] = ode45(@(t,y) ISL1_idealized(t,y,tvec,S4,S2,newopts), tspan, ic);
%         [t,y] = ode45(@(t,y) ISL1_hilln(t,y,tvec,S4,S2,newopts), tspan, ic);
        
        subplot(2,ceil(nconds/2),idx); hold on
        plot(t,y,'LineWidth',3,'Color',[0 0.4470 0.7410 1])
        scatter(tvec(end),ISL1n(cidx),30,[0 0.4470 0.7410],'filled')
        cleanSubplot(12); axis square
        ylim([0,1.1]); xlim([0,tvec(end)])
        legend('calculated','target')
        idx = idx + 1;
        ylabel('ISL1 (au)'); xlabel('time (hrs)')
    end
    sgtitle(sprintf('MSE = %g\n',cmin),'FontSize',14,'FontWeight','bold')
    saveas(gcf,fullfile(savedir,sprintf('ISL1curves_run%.2d.png',ri)))
    
end


%% run the model for all conditions
close all

for ri = 1:nrepeats
paramdir = 'hybrid';
% load(fullfile(savedir,paramdir,sprintf('paramValues_run%.2d',ri)),'newopts')
% load(fullfile(savedir,sprintf('paramValues_run%.2d',ri)),'newopts')
load(fullfile(savedir,sprintf('paramValues_run%.2d',ri)),'thetas','fields','costs')
[cmin,I] = min(costs);
theta = thetas(:,I);
newopts = theta2opts(theta,opts,fields);

% figure('Position',figurePosition(2560,960))
figure('WindowState','maximized')
sISL1 = NaN(ncond,1);
for cidx = 1:size(SOX2,2)
    subplot(3,8,cidx); hold on
    S4 = SMAD4(:,cidx); S2 = SOX2(:,cidx);
    if strcmp(paramdir,'hybrid')
        [t,y] = ode45(@(t,y) ISL1_hybrid(t,y,tvec,S4,S2,newopts), tspan, ic);
    elseif strcmp(paramdir,'ideal')
        [t,y] = ode45(@(t,y) ISL1_idealized(t,y,tvec,S4,S2,newopts), tspan, ic);
    elseif strcmp(paramdir,'doublehill')
        [t,y] = ode45(@(t,y) ISL1_hilln(t,y,tvec,S4,S2,newopts), tspan, ic);
    end
    sISL1(cidx) = y(end);
    
    plot(tvec,S4,'LineWidth',lw)
    plot(tvec,S2,'LineWidth',lw)
    
    plot(t,y,'LineWidth',lw); hold on
    scatter(tvec(end),ISL1n(cidx),100,'filled')
    cleanSubplot(12)
    xlabel('time (hr)'); 
    if mod(cidx - 1,8) == 0
        ylabel('intensity (au)')
    end
    
    if cidx == size(SOX2,2)
        legend('SMAD4','SOX2','sISL1','target','Location','southwest','FontSize',10)
    end
    
    ylim([-0,1.2]); xlim(tvec([1,end]))
end
% saveas(gcf,fullfile(savedir,paramdir,sprintf('ISL1curves_run%.2d.png',ri)))
saveas(gcf,fullfile(savedir,sprintf('allISL1curves_run%.2d.png',ri)))


inds = {[1:8,17],9:17};
colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980];

ISL1s = [ISL1n sISL1];
alphas = [1 0.5]; markers = {'-x','-o'};

figure('Position',figurePosition(560,560)); hold on
for ii = 1:length(inds)
    for jj = 1:size(ISL1s,2)
        if jj == 1
            errorbar(integrals(inds{ii}),ISL1s(inds{ii},jj),ISL1err(inds{ii}),'Color',[colors(ii,:) alphas(jj)],'LineWidth',lw)
        elseif jj == 2
            plot(integrals(inds{ii}),ISL1s(inds{ii},jj),markers{jj},'Color',[colors(ii,:) alphas(jj)],'LineWidth',lw)
        end
    end
end
cleanSubplot(fs); axis square
xlabel('SMAD4 integral'); ylabel('ISL1 (au)')
% saveas(gcf,fullfile(savedir,paramdir,sprintf('ISL1vIntegral_run%.2d.png',ri)))
saveas(gcf,fullfile(savedir,sprintf('ISL1vIntegral_run%.2d.png',ri)))
legend({'42h,target','42h,sim','32h,target','32h,sim'},'FontSize',lfs,'Location','northwest')
end

%% fit SOX2 with inhibition from ISL1
%ISL1 parameters from above fit (should stay static while optimizing SOX2):
% lambdai = 0.124;%0.12;
% alphai = 0.09;%0.11;
% Ksi = 0.266;%0.3315;
% n = 4;

lambdai = 0.106;
alphai = 0.0728;
Ksi = 0.273;
n = 4;

tmax = 42;
ttmax = find(tvec == tmax);
T = tvec(treatmentTime:ttmax)';
ntime = length(T);

tspan = [0,tmax];

%choose a model: linear, hill1, or hill2
model = 'linear';

%take a guess for Kis and ns
Kis = 0.5;
ns = 2;

%initial parameters for regulation of sox2
Ka = 0.85; %autoregulation threshold
na = 4; %autoregulation hill coefficient
betaa = 0.075; %autoregulation amplitude
lambdas = 0.04; %sox2-smad4 slope
alphas = 0.105; %sox2 degradation rate (measured independently)
betas = 0.055; %constitutive production rate

opts = struct('betas',betas,'betaa',betaa,'lambdai',lambdai,'lambda',lambdas,...
    'alphas',alphas,'alphai',alphai,...
    'Ksi',Ksi,'Kis',Kis,'Ka',Ka,...
    'n',n,'ns',ns,'na',na);
fields = {'betas','betaa','lambda','Kis','Ka','ns','na'};

nt = length(fields);

trainconds = 1:17;
nconds = length(trainconds);

%% run the optimization
%create a folder to save results, labeled with the current date
% c = clock;
% savedir = sprintf('%d%.2d%.2d',c(1),c(2),c(3));
% savedir = fullfile(dataDir,'model',savedir(3:end),[dset,'_',s4],'combined');
% if ~exist(savedir,'dir'), mkdir(savedir); end

savedir = fullfile(dataDir,'SOX2_ISL1_autoreg','sigma_em5_temp_em0');
if ~exist(savedir,'dir'), mkdir(savedir); end

%set hyperparameters
mu = zeros(nt,1); %zero vector for mean of gaussian rv
sigma = 5e-5; %variance determines the step size for each parameter 
% Sigma = diag(sigma*[1 1 1 1]); %diagonal covariance matrix
Sigma = diag(sigma*ones(1,nt)); %diagonal covariance matrix
tempmax = 1; %'temperature' determines the acceptance criterion
niters = 1000; %number of iterations per optimization run
nrepeats = 5; %number of independent repeats of the optimzation procedure

for ri = 1:nrepeats
    %randomly initialize parameters
%     theta = rand(nt,1);
    theta = opts2theta(opts,fields);
    newopts = theta2opts(theta,opts,fields);
    V = NaN(ntime,nconds);
    %simulate once to determine initial error
    for idx = 1:nconds
        cidx = trainconds(idx);
        S4 = SMAD4(:,cidx); % %#ok<PFBNS>
%         ic = [mean(SOX2(1:treatmentTime,cidx));0]; % %#ok<PFBNS>
        ic = [1;0];
        %dydt = SOX2_ISL1_hybrid(t,y,tvec,SMAD4,opts)
        % [t,y] = ode45(@(t,y) SOX2_ISL1_hybrid(t,y,tvec,S4,newopts), tspan, ic);
        [t,y] = ode45(@(t,y) SOX2_ISL1_autoreg(t,y,tvec,S4,newopts), tspan, ic);
        V(:,idx) = interp1(t,y(:,1),T);
    end
    
    errold = real(mean(sum((SOX2(treatmentTime:ttmax,trainconds) - V).^2),'all'));
    fprintf('MSE = %g\n',errold)
    
    %initialize variables to store intermediate cost and parameter values
    thetas = zeros(nt,niters);
    costs = zeros(niters,1);
    oldcosts = zeros(niters,1);
    naccept = 0;
    close all
    tic
    %show the error over time
    figure('Position',figurePosition([1250 750]))
    for ii = 1:niters
        %update parameters with random step
        thetap = theta + mvnrnd(mu,Sigma,1)';
        %if any parameters stepped below zero, set them between zero and
        %their initial value instead
        halfp = 0.5*theta;
        thetap(thetap <= 0) = halfp(thetap <= 0);
        newopts = theta2opts(thetap,opts,fields);
        %run the simulation over treatment conditions in parallel
        V = NaN(ntime,nconds);
        parfor idx = 1:nconds
            cidx = trainconds(idx);
            S4 = SMAD4(:,cidx); %#ok<PFBNS>
%             ic = [mean(SOX2(1:treatmentTime,cidx));0]; %#ok<PFBNS>
            ic = [1;0];
            [t,y] = ode45(@(t,y) SOX2_ISL1_autoreg(t,y,tvec,S4,newopts), tspan, ic);
            V(:,idx) = interp1(t,y(:,1),T);
        end
        %error with new parameter values
        errnew = real(mean(sum((SOX2(treatmentTime:ttmax,trainconds) - V).^2),'all'));
        
        costs(ii) = errnew;
        oldcosts(ii) = errold;

        %acceptance criterion
        temp = tempmax*(1 - ii/niters);
        if exp((errold - errnew)/temp) > rand(1)
            theta = thetap;
            errold = errnew;
            naccept = naccept + 1;
        end
        thetas(:,ii) = theta';
        %update the graph of costs every 10 iterations
        if mod(ii,10) == 0
            fprintf('.') %progress indicator
            clf
            semilogy(costs(1:ii),'LineWidth',2)
            hold on
            semilogy(oldcosts(1:ii),'LineWidth',2)
            hold off
            cleanSubplot
            xlabel('iterations'); ylabel('mean squared error')
            axis square
            drawnow
        end
        if mod(ii,500) == 0
            fprintf('\n')
        end
    end
    [cmin,I] = min(costs);
    theta = thetas(:,I);
    newopts = theta2opts(theta,opts,fields);

    fprintf('\nEfficiency = %g\n',naccept/niters)
    fprintf('MSE = %g\n',cmin)
    disp(newopts)
    toc
    
    %for each run, save a graph of proposed and accepted costs over
    %iterations and the associated parameter and cost values as mat files
    close all
    figure('Position',figurePosition([1250 750]))
    semilogy(costs,'LineWidth',2); hold on
    semilogy(oldcosts,'LineWidth',2); hold off
    xlabel('iterations'); ylabel('mean squared error')
    cleanSubplot; axis square
%     ylim([0,4])
    legend('proposed cost', 'accepted cost')
    saveas(gcf,fullfile(savedir,sprintf('costs_run%.2d.png',ri)))
    
    
    % newopts = theta2opts(theta,opts,fields);
    save(fullfile(savedir,sprintf('paramValues_run%.2d',ri)),'newopts','thetas','fields','costs','oldcosts')
    
    margin = 0.02;
    close all
    figure('Position',figurePosition(2500, 700))
    idx = 1;
    for ii = 1:nconds
        
        cidx = trainconds(ii);
        ytest = SOX2(treatmentTime:ttmax,cidx); S4 = SMAD4(:,cidx);
        
%         ic = [mean(SOX2(1:treatmentTime,cidx));0];
        ic = [1;0];
        [t,y] = ode45(@(t,y) SOX2_ISL1_autoreg(t,y,tvec,S4,newopts), tspan, ic);
        y = real(y);
        
        subplot(2,ceil(nconds/2),idx); hold on
        %SOX2
        plot(t,y(:,1),'LineWidth',3,'Color',[0 0.4470 0.7410 0.5])
        plot(T,ytest,'LineWidth',1.5,'Color',[0 0.4470 0.7410 1])
        %ISL1
        plot(t,y(:,2),'LineWidth',1.5,'Color',[0.8500 0.3250 0.0980 0.5])
        scatter(tvec(end),ISL1n(cidx),30,[0.8500 0.3250 0.0980],'filled')
        
        cleanSubplot(12); axis square
        ylim([0,1.1]); xlim([0,tvec(end)])
        if ii == nconds
            legend('calculated SOX2','target SOX2','calculated ISL1','target ISL1',...
                'Location','southeast')
        end
        
        idx = idx + 1;
        ylabel('SOX2::GFP'); xlabel('time (hrs)')
    end
    sgtitle(sprintf('MSE = %g\n',cmin),'FontSize',14,'FontWeight','bold')
    saveas(gcf,fullfile(savedir,sprintf('SOX2_ISL1_curves_run%.2d.png',ri)))
    
end
% t = NaN; y = NaN; ic = NaN;

% compare model to experiment for all conditions
for ri = 1:nrepeats
load(fullfile(savedir,sprintf('paramValues_run%.2d',ri)),'newopts')

figure('Position',figurePosition(2560,960))
for ii = 1:size(SOX2,2)
        
    cidx = ii;
    ytest = SOX2(:,cidx); S4 = SMAD4(:,cidx);
        
    ic = [mean(SOX2(1:treatmentTime,cidx));0];
%     ic = [1;0];
    [t,y] = ode45(@(t,y) SOX2_ISL1_autoreg(t,y,tvec,S4,newopts), tspan, ic);
    y = real(y);

    subplot(3,8,ii); hold on
    %SOX2
    plot(t,y(:,1),'LineWidth',3,'Color',[0 0.4470 0.7410 0.5])
    plot(tvec,ytest,'LineWidth',1.5,'Color',[0 0.4470 0.7410 1])
    %ISL1
    plot(t,y(:,2),'LineWidth',1.5,'Color',[0.8500 0.3250 0.0980 0.5])
    scatter(tvec(end),ISL1n(cidx),30,[0.8500 0.3250 0.0980],'filled')

    cleanSubplot(12); axis square
    ylim([0,1.1]); xlim([0,tvec(end)])
    if ii == ncond
        legend('calculated SOX2','target SOX2','calculated ISL1','target ISL1',...
            'Location','southeast')
    end
    title(meta.conditions(ii))
    ylabel('SOX2::GFP'); xlabel('time (hrs)')
end
saveas(gcf,fullfile(savedir,sprintf('all_SOX2_ISL1_curves_run%.2d.png',ri)))
end

%% plot simulated and measured SOX2 vs time by duration
% close all
% paramdir = 'D:\220905_SOX2_Smad4GFP_BMP_IWP2_LDN_live\model\230206\220830_idealS4\combined';
% paramname = 'paramValues_run01.mat';
% paramdir = 'D:\220905_SOX2_Smad4GFP_BMP_IWP2_LDN_live\model\230214\220830_idealS4\combined';
% paramname = 'paramValues_run01.mat';

% paramdir = 'G:\My Drive\Research\Heemskerk lab\Data\220830_220905_SOX2GFP_model\avg_SOX2_idealS4_SMAD4\combined';
% paramname = 'paramValues_run04.mat';

paramdir = fullfile(dataDir,'SOX2_ISL1_autoreg','sigma_em5_temp_em0');
paramname = 'paramValues_run05.mat';

% load(fullfile(paramdir,paramname),'newopts')
load(fullfile(paramdir,paramname),'thetas','fields','costs')
[~,I] = min(costs);
theta = thetas(:,I);
newopts = theta2opts(theta,opts,fields);
tspan = [0 42];

conds = {[1:8,17],[9:16,17]};
% ax = gobjects(length(conds),1);
legpos = [0.725 0.3634 0.2286 0.5563];
legstr = meta.conditions(conds{1});
legstr = strrep(legstr,'LDN',''); legstr = strrep(legstr,',42hr',''); legstr = strrep(legstr,'mTeSR','X');
titles = {'42 hours','32 hours'}; savelabels = {'42hr','32hr'};

sISL1 = NaN(length(conds{1}),length(conds));
mISL1 = NaN(length(conds{1}),length(conds));
errISL1 = NaN(length(conds{1}),length(conds));
for ii = 1:length(conds)
    colors = turbo(length(conds{ii}));
    colors = colors(end:-1:1,:);
    figure('Position',figurePosition([560 560]))
    hold on
    p = gobjects(length(conds{ii}),1);
    for jj = 1:length(conds{ii})
        cidx = conds{ii}(jj);
        S4 = SMAD4(:,cidx);
        % ic = [mean(SOX2(1:treatmentTime,cidx));0];
        ic = [1;0];
        [t,y] = ode45(@(t,y) SOX2_ISL1_autoreg(t,y,tvec,S4,newopts), tspan, ic);
        plot(t,y(:,1),'LineWidth',3,'Color',[colors(jj,:),0.5])
        p(jj) = plot(tvec,SOX2(:,cidx),'LineWidth',1.5,'Color',[colors(jj,:),1]);
        sISL1(jj,ii) = y(end,2);
        mISL1(jj,ii) = ISL1n(cidx);
        errISL1(jj,ii) = ISL1err(cidx);
    end
    hold off
    cleanSubplot(fs); axis square
    xlim([min(tvec),max(tvec)])
    
    ylabel('GFP::SOX2 (au)')
    xlabel('time (hr)')
%     ylim([-0.05,1.05])%ylim([400,3300])%ylim([600,1800])
    ylim([min(SOX2,[],'all'), max(SOX2,[],'all')] + [-0.05,0.05])
    
    savefigure(fullfile(dataDir,'figures',...
        strcat('SOX2_',savelabels{ii},'_datavmodel.png')))
    savefigure(fullfile(paramdir,...
        strcat('SOX2_',savelabels{ii},'_datavmodel.png')))
    
    if ii == 1
        lgd = legend(p,legstr,'Position',legpos,'FontSize',lfs);
        title(lgd,{'BMPRi','(nM)'},'FontSize',lfs)
        savefigure(fullfile(dataDir,'figures',...
            strcat('SOX2_',savelabels{ii},'_datavmodel_withlegend.png')))
        savefigure(fullfile(paramdir,...
            strcat('SOX2_',savelabels{ii},'_datavmodel_withlegend.png')))
    end
end

ISL1s = {sISL1,mISL1};
ylabels = {'simulated ISL1 (au)','ISL1 (au)'};
legstr = {'42hr','32hr'}; 
savelabels = {'ISL1simulated','ISL1measured'};
%add error bars to measured ISL1?
for jj = 1:length(ISL1s)
    figure('Position',figurePosition(560,560)); hold on
    for ii = 1:length(conds)
        if jj == 1
            plot(integrals(conds{ii}),ISL1s{jj}(:,ii),'-x','LineWidth',lw)
        elseif jj == 2
            errorbar(integrals(conds{ii}),ISL1s{jj}(:,ii),ISL1err(conds{ii}),'LineWidth',lw)
        end
    end
    legend(legstr,'Location','northwest','FontSize',lfs)
    cleanSubplot(fs); axis square
    xlabel('SMAD4 (N:C) integral')
    ylabel(ylabels{jj})
    
    savefigure(fullfile(dataDir,'figures',...
        strcat(savelabels{jj},'_v_integral.png')))
    savefigure(fullfile(paramdir,...
        strcat(savelabels{jj},'_v_integral.png')))
end

%% local functions

function opts = theta2opts(theta,opts,fields)
if length(theta) ~= length(fields)
    error('mismatch in number of parameters')
end

for ii = 1:length(theta)
    opts.(fields{ii}) = theta(ii);
end

end


function theta = opts2theta(opts,fields)

theta = NaN(length(fields),1);
for ii = 1:length(fields)
    theta(ii) = opts.(fields{ii});
end


end

function X = preprocess_SMAD4_220629(S4struct)

signalData = S4struct.signalData;
tvec = S4struct.tvec;
ntime = length(tvec);

maxcond = 1;

X = cellfun(@(x) mean(x,'omitnan'),signalData);
%normalize
maxval = mean(X(:,maxcond));
minval = mean(X(tvec <= 0,:),'all');
X = (X - minval)/(maxval - minval);

x1 = X(:,1:8);
x2 = X(:,9:16);

%interpolate to infer signaling in the 100 nM LDN dose condition
ds = log10(1 + [0 1 5 10 25 50 250 1000]);
dq = log10(1 + 100);
v1 = NaN(ntime,1);
for ti = 1:ntime
    vs = x1(ti,:);
    vq = interp1(ds,vs,dq);
    v1(ti) = vq;
end

ds = log10(1 + [0 1 5 10 25 50 250 1000]);
dq = log10(1 + 100);
v2 = NaN(ntime,1);
for ti = 1:ntime
    vs = x2(ti,:);
    vq = interp1(ds,vs,dq);
    v2(ti) = vq;
end

X = [X(:,[1,3:6]),v1,X(:,7:8),X(:,[9,11:14]),v2,X(:,15:16),zeros(ntime,1)];

end


function SMAD4 = makeIdealS4(tvec,s4levels)

ntime = length(tvec);
nlevel = length(s4levels)-1;
l2 = s4levels(nlevel);

lambda = 2;

SMAD4 = zeros(ntime,2*nlevel+1);
for ii = 1:nlevel
    t1 = tvec >= 0;
    t2 = tvec >= 32;
    l1 = s4levels(ii);
    
    v1 = s4levels(ii)*(1 - exp(-lambda*tvec(t1)));
    v2 = l2 + (l1 - l2)*exp(-lambda*(tvec(t2) - 32));
    
    SMAD4(t1, ii) = v1;
    SMAD4(t1, nlevel + ii) = v1;
    SMAD4(t2, nlevel + ii) = v2;
end
end



















