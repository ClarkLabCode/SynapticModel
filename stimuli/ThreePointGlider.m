function [ stimArray ] = ThreePointGlider(param, noiseParam, updateRate)
%UNTITLED16 Summary of this function goes here
%   Detailed explanation goes here

% Compute the number of discrete points in space and time
numX = floor(param.xTot / noiseParam.barWidth);
numT = floor(param.tTot * updateRate);

% Make the 1-dimensional seed image
seed = 2 * double(rand(1, numX) > 0.5) - 1;

% Allocate containers
PDcon = nan(numT,numX);
NDcon = nan(numT,numX);

% Initialize each container
PDcon(1,:) = seed;
NDcon(1,:) = seed;

% Fill in each kymograph timestep by timestep
for ind = 2:numT
    PDcon(ind,:) = circshift(PDcon(ind-1,:), +1, 2) .* circshift(PDcon(ind-1,:), +0, 2);
    NDcon(ind,:) = circshift(NDcon(ind-1,:), -1, 2) .* circshift(NDcon(ind-1,:), +0, 2);
end

% Make the diverging gliders the easy way
PDdiv = flipud(NDcon);
NDdiv = flipud(PDcon);

% Make the uncorrelated control
UN = 2 * double(rand(numT, numX) > 0.5) - 1;

% Resample to the appropriate size using nearest-neighbor interpolation
PDcon = imresize(PDcon, [length(param.t), length(param.x)], 'nearest');
NDcon = imresize(NDcon, [length(param.t), length(param.x)], 'nearest');
PDdiv = imresize(PDdiv, [length(param.t), length(param.x)], 'nearest');
NDdiv = imresize(NDdiv, [length(param.t), length(param.x)], 'nearest');
UN    = imresize(UN, [length(param.t), length(param.x)], 'nearest');

% Combine everything together
stimArray = noiseParam.mlum + noiseParam.c * param.mask .* cat(3, PDcon, -PDcon, NDcon, -NDcon, PDdiv, -PDdiv, NDdiv, -NDdiv, UN);

end
