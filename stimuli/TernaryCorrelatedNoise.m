function [ stimArray ] = TernaryCorrelatedNoise(param, noiseParam, dt)
%TERNARYCORRELATEDNOISE: Ternary noise with pairwise correlations.

% Compute the number of discrete points in space and time
numX = floor(param.xTot / noiseParam.barWidth);
numT = floor(param.tTot * noiseParam.updateRate);

% Make the binary noise arrays
B = 2 * double(rand(numT, numX) > 0.5) - 1;
C = 2 * double(rand(numT, numX) > 0.5) - 1;

% Make the stimuli
PosCorr = (B + circshift(B, -[dt, noiseParam.dx]))/2;
NegCorr = (B - circshift(B, -[dt, noiseParam.dx]))/2;
UnCorr  = (B + C)/2;

% Resample to the appropriate size using nearest-neighbor interpolation
PosCorr = imresize(PosCorr, [length(param.t), length(param.x)], 'nearest');
NegCorr = imresize(NegCorr, [length(param.t), length(param.x)], 'nearest');
UnCorr  = imresize(UnCorr,  [length(param.t), length(param.x)], 'nearest');

% Combine everything together
stimArray = noiseParam.mlum + noiseParam.c * param.mask .* cat(3, PosCorr, NegCorr, UnCorr);

end

