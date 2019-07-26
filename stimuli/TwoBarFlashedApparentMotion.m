function [ stimArray ] = TwoBarFlashedApparentMotion(params, barParam, barCenter)

% Grab parameters
w = barParam.barWidth;
d = barParam.duration;
t = params.t;
x = params.x;

% Allocate a container
seed = zeros(length(t),length(x),16);

% ON components of PD sequence
seed(:,:,1) = makeFlashedBar(t, x, d, 0*d, w, barCenter - w);
seed(:,:,2) = makeFlashedBar(t, x, d, 1*d, w, barCenter + 0);

% ON components of ND sequence
seed(:,:,3) = makeFlashedBar(t, x, d, 0*d, w, barCenter + w);
seed(:,:,4) = makeFlashedBar(t, x, d, 1*d, w, barCenter + 0);

% OFF components
seed(:,:,5:8) = -seed(:,:,1:4);

% Sequences: {'++PD','++ND','--PD','--ND','+-PD','+-ND','-+PD','-+ND'}
seed(:,:,9)  = sum(seed(:,:,[1,2]),3);
seed(:,:,10) = sum(seed(:,:,[3,4]),3);
seed(:,:,11) = sum(seed(:,:,[5,6]),3);
seed(:,:,12) = sum(seed(:,:,[7,8]),3);
seed(:,:,13) = sum(seed(:,:,[1,6]),3);
seed(:,:,14) = sum(seed(:,:,[3,8]),3);
seed(:,:,15) = sum(seed(:,:,[2,5]),3);
seed(:,:,16) = sum(seed(:,:,[4,7]),3);

% Adjust luminance and contrast
stimArray = barParam.mlum + barParam.c * params.mask .* seed;

end

function [ seed ] = makeFlashedBar(t, x, dur, tau, barWidth, barPos)
seed = double((t >= tau) & (t < tau + dur) & ((x >= barPos-barWidth/2) & (x < barPos+barWidth/2)));
end