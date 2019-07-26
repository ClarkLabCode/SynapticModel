function ThreeInputModelSinusoidTimetraces()

%% Set parameters

% Set overall parameters
[ params ] = SetModelParameters('tInt', 0, 'tOn', 4, 'dx', 1, 'tAv', 2);

% Make the filters
[ filters ] = MakeModelFilters(params);

% Contrast
c = 1/2;

% Temporal and spatial frequency vectors
tf = 1;
sf = 1/45;

%% Compute everything

% Find spatial central point
[~,spatialInd] = min(abs(params.x - 180));

% Generate gratings
[ stimArray ] = MonocomponentSinusoidalGratings(params, c, tf, sf);

% Compute photoreceptor responses
[ resp1, resp2, resp3 ] = ComputePhotoreceptorResponses(stimArray, params, filters);

% Define the input rectifiers
if isinf(params.inputRectBeta)
    relu = @(f) (f .* (f>0));
else
    relu = @(f) f.*(erf(params.inputRectBeta * f)+1)/2;
end

% Define the output half-quadratic
if isinf(params.outputRectBeta)
    halfsquare = @(f) (f .* (f>0)).^2;
else
    halfsquare = @(f) (f.*(erf(params.outputRectBeta*f)+1)/2).^2;
end

% Compute each postsynaptic conductance
% Input 1 ~ Mi9
% Input 2 ~ Mi1
% Input 3 ~ Mi4
g1 = params.g1 .* relu(-resp1);
g2 = params.g2 .* relu(+resp2);
g3 = params.g3 .* relu(+resp3);

% Compute the numerator and denominator of the three-input model
numResp = params.V1 .* g1 + params.V2 .* g2 + params.V3 .* g3;
denResp = params.gleak + g1 + g2 + g3;

% Compute the voltage response of the full model
voltageResp = numResp ./ denResp;

% Model the transformation from membrane voltage to calcium concentration
% as a half-quadratic
calciumResp = halfsquare(voltageResp);

%% Plot the results

t = params.t(params.averagingMask) - params.tAv;

MakeFigure;
subplot(4,1,1);
prShift = floor(params.photoreceptorSpacing / params.dx);
respMat = stimArray(params.averagingMask,spatialInd + (-1:1) * prShift,1);
plot(t, respMat,'linewidth', 2);
ylabel('contrast');
xlabel('time (s)');

subplot(4,1,2);
respMat = [g1(params.averagingMask,spatialInd,1), g2(params.averagingMask,spatialInd,1), g3(params.averagingMask,spatialInd,1)];
plot(t, respMat,'linewidth', 2);
ylabel('conductance/g_leak');

subplot(4,1,3);
respMat = voltageResp(params.averagingMask,spatialInd,1);
plot(t, respMat,'linewidth', 2);
ylabel('V_{mem} (mV)');

subplot(4,1,4);
respMat = calciumResp(params.averagingMask,spatialInd,1);
plot(t, respMat,'linewidth', 2);
ylabel('response (arb. units)');

end