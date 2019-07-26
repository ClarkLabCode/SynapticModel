function [ stimArray ] = FullFieldBinaryFlicker(params, noiseParam, c)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% Compute the number of discrete points in time
numT = floor(params.tTot * noiseParam.updateRate);

% Make the 1-dimensional seed image
seed = 2 * double(rand(numT, 1) > 0.5) - 1;

% Resample to the appropriate size using nearest-neighbor interpolation
seed = imresize(seed, [length(params.t), 1], 'nearest');

% Expand to fill the full space
seed = repmat(seed, 1, length(params.x));

% Adjust luminance and contrast
stimArray = noiseParam.mlum + c * params.mask .* seed;

end

