function [ stimArray ] = CompositeSinusoidalGratingsForLinearityAnalysis(params, c, tf, sf)
% Composite sinusoidal gratings for linearity analysis
% Output format is [time x space x numStim x phaseShift]

% Convert time and space vectors into units of radians
tVec = 2 * pi * tf * params.t;
xVec = 2 * pi * sf * params.x;

% Phase offsets used in Wienecke et al. 2018
phi = (0:1:7)/8*pi; 

% Reshape phase shift vectors appropriately
phi = permute(phi, [4,3,1,2]);

% Make the stimuli

% Allocate a container
stimArray = zeros(length(tVec), length(xVec), 2, length(phi));

% PD composite
stimArray(:,:,1,:) = sin(tVec + phi - pi/2) .* sin(xVec + phi);

% ND composite
stimArray(:,:,2,:) = sin(tVec + phi + pi/2) .* sin(xVec - phi);

% Adjust contrast and mask stimulus
stimArray = c .* params.mask .* stimArray;

end

