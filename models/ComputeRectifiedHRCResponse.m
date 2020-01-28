function [meanResp, voltageResp, calciumResp ] = ComputeRectifiedHRCResponse(stimArray, p, f)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%% Compute photoreceptor responses

% Blur the filter in space (assumes periodic boundary conditions)
if p.useSpatialFilter
    blurArray = fftshift(ifft(fft(f.spatialFilter,[],2) .* fft(stimArray,[],2), [], 2),2);
else
    blurArray = stimArray;
end

% Filter the spatially-blurred stimulus in time
lpArray = filter(f.lp, 1, blurArray, [], 1);
hpArray = filter(f.hp, 1, blurArray, [], 1);

% Shift the stimulus to get each of the photoreceptor inputs
% Note that signs of circshifts are reversed relative to index notation
prShift = floor(p.photoreceptorSpacing / p.dx);
resp1 = circshift(lpArray, +prShift, 2);
resp2 = hpArray;
resp3 = circshift(hpArray, +prShift, 2);
resp4 = lpArray;


%% Compute HRC response

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

% Standard HRC
voltageResp = resp1 .* resp2 - resp3 .* resp4;

% Rectify output
calciumResp = relu(voltageResp);
% calciumResp = voltageResp;

%% Average model responses over time and phase

meanResp = squeeze(nanmean(nanmean(nanmean(calciumResp(p.averagingMask,:,:),4),2),1));

end

