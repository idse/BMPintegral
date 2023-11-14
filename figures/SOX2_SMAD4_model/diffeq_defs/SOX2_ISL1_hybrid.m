function dydt = SOX2_ISL1_hybrid(t,y,tvec,SMAD4,opts)
%t = current time
%y = level of each transciption factor at time t
%components of y:
%y(1) = SOX2; y(2) = ISL1
%tvec = vector of time points at which SMAD4 levels are sampled
%SMAD4 = SMAD4 levels at each time in tvec
%opts = structure array of parameters for the ODE model

%interpolate BMP and Activin levels at the given time point
if numel(tvec) > 1
    smad4 = interp1(tvec,SMAD4,t,'nearest'); % Interpolate at time t
elseif numel(tvec) == 1 %if only one level is provided, no need to interpolate
    smad4 = SMAD4;
end

%%%parameters%%%
%production rates
betas = opts.betas; lambdai = opts.lambdai;
%degradation rates
alphas = opts.alphas; alphai = opts.alphai;
%hill function thresholds
Ksi = opts.Ksi; Kis = opts.Kis;
%SMAD4/SOX2 slope
lambda = opts.lambda;
%hill coefficient for regulation of ISL1 by SOX2
n = opts.n;
%hill coefficient for regulation of SOX2 by ISL1
ns = opts.ns;

dydt = NaN(size(y));
%sox2:
dydt(1) = (betas - lambda*smad4)/(1 + (y(2)/Kis)^ns) - alphas*y(1);
%isl1:
dydt(2) = lambdai*smad4/(1 + (y(1)/Ksi)^n) - alphai*y(2);


end