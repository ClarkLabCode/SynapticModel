function [ stimArray ] = UncorrelatedBinaryNoise(param, noiseParam)

% Compute the number of discrete points in space and time
numX = floor(param.xTot / noiseParam.barWidth);
numT = floor(param.tTot * noiseParam.updateRate);

% Make the binary noise array
B = 2 * double(rand(numT, numX) > 0.5) - 1;

% Resample to the appropriate size using nearest-neighbor interpolation
noiseArray = imresize(B, [length(param.t), length(param.x)], 'nearest');

% Combine everything together
stimArray = noiseParam.mlum + noiseParam.c * param.mask .* noiseArray;

end

