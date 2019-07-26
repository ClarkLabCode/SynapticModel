function [ stimArray ] = OrientedSinusoidGratings(params, filters, c, tf, sf, theta)
%Oriented sinusoidal gratings, filtered in second spatial dimension.
% Assumes that input theta is in degrees.

% % Reformat spacetime vectors to dimensions [y,t,x]
% y = permute(params.x,[1,3,2]);
% k = permute(filters.spatialFilter, [1,3,2]);
% 
% stimArray = nan(length(params.t),length(params.x),length(theta));
% 
% for indT = 1:length(theta)
%     
%     % Compute time and space vectors
%     xVec = 2 * pi * sf * (cos(theta(indT)).* params.x + sin(theta(indT)) .* y);
%     tVec = 2 * pi * tf * params.t;
%     
%     % Compute the stimulus and contract using the spatial filter to dimensions [t,x]
%     stimArray(:,:,indT) = squeeze(sum(k .* sin(tVec - xVec),3));
% end
% 
% % Adjust contrast and mask stimulus
% stimArray = c .* params.mask .* stimArray;


% Compute time and space vectors
tVec = 2 * pi * tf * params.t;
xVec = 2 * pi * sf * cos(theta) * params.x;

% Compute amplitude
sBlur = params.fwhmBlur / (2*sqrt(2*log(2)));
a = exp(-2*(pi * sf * sBlur * sin(theta)).^2);

% Adjust contrast and mask stimulus
stimArray = c .* a .* params.mask .* sin(tVec-xVec);

end

