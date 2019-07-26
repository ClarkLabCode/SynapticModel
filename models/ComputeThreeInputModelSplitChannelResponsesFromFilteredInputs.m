function [calciumResp] = ComputeThreeInputModelSplitChannelResponsesFromFilteredInputs(resp1, resp2, resp3, p)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% Define the input rectifiers
if isinf(p.inputRectBeta)
    relu = @(f) (f .* (f>0));
else
    relu = @(f) f.*(erf(p.inputRectBeta * f)+1)/2;
end

% Define the output half-quadratic
if isinf(p.outputRectBeta)
    halfsquare = @(f) (f .* (f>0)).^2;
else
    halfsquare = @(f) (f.*(erf(p.outputRectBeta*f)+1)/2).^2;
end

%% T4 progressive

% Compute each postsynaptic conductance
g1 = p.g1 .* relu(-resp1);
g2 = p.g2 .* relu(+resp2);
g3 = p.g3 .* relu(+resp3);

% Compute the numerator and denominator of the three-input model
numResp = p.V1 .* g1 + p.V2 .* g2 + p.V3 .* g3;
denResp = p.gleak + g1 + g2 + g3;

% Compute the voltage response of the full model
voltageResp = numResp ./ denResp;

% Model the transformation from membrane voltage to calcium concentration
% as a half-quadratic
t4ProCalciumResp = halfsquare(voltageResp);

%% T4 regressive

% Compute each postsynaptic conductance
g1 = p.g1 .* relu(-resp3);
g2 = p.g2 .* relu(+resp2);
g3 = p.g3 .* relu(+resp1);

% Compute the numerator and denominator of the three-input model
numResp = p.V1 .* g1 + p.V2 .* g2 + p.V3 .* g3;
denResp = p.gleak + g1 + g2 + g3;

% Compute the voltage response of the full model
voltageResp = numResp ./ denResp;

% Model the transformation from membrane voltage to calcium concentration
% as a half-quadratic
t4RegCalciumResp = halfsquare(voltageResp);

%% T5 progressive

% Compute each postsynaptic conductance
g1 = p.g1 .* relu(+resp1);
g2 = p.g2 .* relu(-resp2);
g3 = p.g3 .* relu(-resp3);

% Compute the numerator and denominator of the three-input model
numResp = p.V1 .* g1 + p.V2 .* g2 + p.V3 .* g3;
denResp = p.gleak + g1 + g2 + g3;

% Compute the voltage response of the full model
voltageResp = numResp ./ denResp;

% Model the transformation from membrane voltage to calcium concentration
% as a half-quadratic
t5ProCalciumResp = halfsquare(voltageResp);

%% T5 regressive

% Compute each postsynaptic conductance
g1 = p.g1 .* relu(+resp3);
g2 = p.g2 .* relu(-resp2);
g3 = p.g3 .* relu(-resp1);

% Compute the numerator and denominator of the three-input model
numResp = p.V1 .* g1 + p.V2 .* g2 + p.V3 .* g3;
denResp = p.gleak + g1 + g2 + g3;

% Compute the voltage response of the full model
voltageResp = numResp ./ denResp;

% Model the transformation from membrane voltage to calcium concentration
% as a half-quadratic
t5RegCalciumResp = halfsquare(voltageResp);

%% Combine all responses together

calciumResp = cat(5, t4ProCalciumResp, t4RegCalciumResp, t5ProCalciumResp, t5RegCalciumResp);

end

