function [ stimArray ] = MonocomponentSinusoidalGratings(params, c, tf, sf)
% Simple monocomponent sinusoidal gratings
% Output format is [time x space x numStim]

% Convert time and space vectors into units of radians
tVec = 2 * pi * tf * params.t;
xVec = 2 * pi * sf * params.x;

% Allocate a container
stimArray = zeros(length(tVec), length(xVec), 2);

% PD
stimArray(:,:,1) = sin(tVec - xVec);

% ND
stimArray(:,:,2) = sin(tVec + xVec);

% Adjust contrast and mask stimulus
stimArray = c .* params.mask .* stimArray;

end