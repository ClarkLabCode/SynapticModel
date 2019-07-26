function [ stimArray ] = TwoPointGlider(param, noiseParam, updateRate)
%UNTITLED16 Summary of this function goes here
%   Detailed explanation goes here

% Compute the number of discrete points in space and time
numX = floor(param.xTot / noiseParam.barWidth);
numT = floor(param.tTot * updateRate);

% Make the 1-dimensional seed image
seed = 2 * double(rand(1, numX) > 0.5) - 1;

% Allocate containers
PDpos = nan(numT,numX);
NDpos = nan(numT,numX);
PDneg = nan(numT,numX);
NDneg = nan(numT,numX);
% UN = nan(numT,numX);

% Initialize each container
PDpos(1,:) = seed;
NDpos(1,:) = seed;
PDneg(1,:) = seed;
NDneg(1,:) = seed;

% Fill in each kymograph timestep by timestep
for ind = 2:numT
    PDpos(ind,:) = +circshift(PDpos(ind-1,:), +1, 2);
    NDpos(ind,:) = +circshift(NDpos(ind-1,:), -1, 2);
    PDneg(ind,:) = -circshift(PDneg(ind-1,:), +1, 2);
    NDneg(ind,:) = -circshift(NDneg(ind-1,:), -1, 2);
end

% Make the uncorrelated control
UN = 2 * double(rand(numT, numX) > 0.5) - 1;

% Resample to the appropriate size using nearest-neighbor interpolation
PDpos = imresize(PDpos, [length(param.t), length(param.x)], 'nearest');
NDpos = imresize(NDpos, [length(param.t), length(param.x)], 'nearest');
PDneg = imresize(PDneg, [length(param.t), length(param.x)], 'nearest');
NDneg = imresize(NDneg, [length(param.t), length(param.x)], 'nearest');
UN    = imresize(UN,    [length(param.t), length(param.x)], 'nearest');

% Combine everything together
stimArray = noiseParam.mlum + noiseParam.c * param.mask .* cat(3, PDpos, NDpos, PDneg, NDneg, UN);

end
