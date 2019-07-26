function [ stimArray ] = LocalFlashedBar(params, barParam, dur)

% Permute durations
dur = permute(dur, [3,2,1]);

% Make the seed image
locX = mean(params.x);
seed = (params.t >= 0) & (params.t <= dur) & ((params.x >= locX-barParam.barWidth/2) & (params.x <= locX+barParam.barWidth/2));
seed = double(seed);

% Adjust luminance and contrast
stimArray = barParam.mlum + barParam.c * params.mask .* cat(3, seed,-seed);

end