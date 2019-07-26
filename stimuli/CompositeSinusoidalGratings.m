function [ stimArray ] = CompositeSinusoidalGratings(params, c, tf, sf)
% PD + ND and PD + OD composite sinusoidal gratings
% Output format is [time x space x numStim x phaseShift]

% Convert time and space vectors into units of radians
tVec = 2 * pi * tf * params.t;
xVec = 2 * pi * sf * params.x;

% Define phase offsets
numShift = params.numRep;
if params.useRandomShifts
    % Sample phase shifts randomly
    phi1 = 2 * pi * [0, rand(1, numShift-1)];
    phi2 = 2 * pi * [0, rand(1, numShift-1)];
else    
    % Sample numShift points from the grid of phase shifts
    nn = floor(sqrt(numShift));
    numShift = nn^2;
    [phi1, phi2] = meshgrid((0:1:nn-1)*2*pi/nn);
    phi1 = phi1(:)';
    phi2 = phi2(:)';
end

% Reshape phase shift vectors appropriately
phi1 = permute(phi1, [4,3,1,2]);
phi2 = permute(phi2, [4,3,1,2]);

%% Make the stimuli

% Allocate a container
stimArray = zeros(length(tVec), length(xVec), 2, numShift);

% PD + ND
stimArray(:,:,1,:) = sin(tVec - (xVec + phi1)) + sin(tVec + (xVec + phi2));

% PD + OD
stimArray(:,:,2,:) = sin(tVec - (xVec + phi1)) + sin(tVec + phi2);

% Adjust contrast and mask stimulus
stimArray = c .* params.mask .* stimArray;

end

