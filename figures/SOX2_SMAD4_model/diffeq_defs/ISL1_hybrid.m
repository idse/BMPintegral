function dydt = ISL1_hybrid(t,y,tvec,SMAD4,SOX2,opts)
%t = current time
%y = level of ISL1 at time t
%tvec = vector of time points at which SMAD4 and SOX2 levels are sampled
%SMAD4 = SMAD4 levels at each time in tvec
%SOX2 = SOX2 levels at each time in tvec
%opts = structure array of parameters for the ODE model
%a high-order hill function governs regulation of ISL1 by SOX2

%interpolate BMP and Activin levels at the given time point
if numel(tvec) > 1
    smad4 = interp1(tvec,SMAD4,t,'nearest'); % Interpolate at time t
    sox2 = interp1(tvec,SOX2,t,'nearest');
elseif numel(tvec) == 1 %if only one level is provided, no need to interpolate
    smad4 = SMAD4;
    sox2 = SOX2;
end

%%%parameters%%%
%production rate
% betai = opts.betai;
%degradation rate
alphai = opts.alphai;
%hill function thresholds
Ksi = opts.Ksi; lambdai = opts.lambdai;
%hill coefficient (change to different values for different functions?)
n = opts.n;

dydt = lambdai*smad4/(1 + (sox2/Ksi)^n) - alphai*y;

end