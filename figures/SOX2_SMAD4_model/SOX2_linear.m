function dydt = SOX2_linear(t,y,tvec,SMAD4,opts)
%t = current time
%y = level SOX2 at time t
%tvec = vector of time points at which SMAD4 levels are sampled
%SMAD4 = SMAD4 levels at each time in tvec
%opts = structure array of parameters for the ODE model

%interpolate SMAD4 level at the given time point
if numel(tvec) > 1
    smad4 = interp1(tvec,SMAD4,t,'nearest'); % Interpolate at time t
elseif numel(tvec) == 1 %if only one level is provided, no need to interpolate
    smad4 = SMAD4;
end

%%%parameters%%%
%production rate
beta = opts.beta;
%degradation rate
alpha = opts.alpha;
%SMAD4/SOX2 slope
lambda = opts.lambda;

dydt = beta - lambda*smad4 - alpha*y;


end