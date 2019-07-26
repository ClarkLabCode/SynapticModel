function [ stimArray ] = RandomCheckerboard(p, noiseParam, v)

% Set the number of stimuli and extract needed parameters
numStim = 6;
rho = noiseParam.density;
numX = floor(p.xExtent / noiseParam.barWidth);

% Make the 1-dimensional seed image
eta = rand(numX, numX);
seed = imresize(double(eta <= rho) - double(eta > (1-rho)),[length(p.x), length(p.x)], 'nearest');

% Allocate a container
stimArray = zeros(length(p.t), length(p.x), numStim);

% Compute the offsets at each timestep
dx = floor(p.t * v * p.dx);

% Randomize sign of OD motion
odSign = sign(2*rand(1) - 1);

% Shift the seed the appropriate amount
for ind = 1:length(p.t)
    
    % PD
    temp = circshift(seed, +dx(ind), 2);
    stimArray(ind,:,1) = temp(1,:);
    
    % ND
    temp = circshift(seed, -dx(ind), 2);
    stimArray(ind,:,2) = temp(1,:);
    
    % OD
    temp = circshift(seed, odSign*dx(ind), 1);
    stimArray(ind,:,6) = temp(1,:);
    
    % PD + ND
    stimArray(ind,:,3) = clippedsum(stimArray(ind,:,1), stimArray(ind,:,2));
    
    % PD + OD
    stimArray(ind,:,4) = clippedsum(stimArray(ind,:,1), stimArray(ind,:,6));
    
    % ND + OD
    stimArray(ind,:,5) = clippedsum(stimArray(ind,:,2), stimArray(ind,:,6));
end

% Mask and adjust the stimulus
stimArray = noiseParam.mlum + noiseParam.c * p.mask .* stimArray;

end


function [ z ] = clippedsum(x, y)
z = double((x>0) | (y>0)) - double((x<0) | (y<0));
end
