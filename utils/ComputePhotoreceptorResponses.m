function [ resp1, resp2, resp3, blurArray ] = ComputePhotoreceptorResponses(stimArray, p, f)

% Blur the filter in space (assumes periodic boundary conditions)
if p.useSpatialFilter
    blurArray = fftshift(ifft(fft(f.spatialFilter,[],2) .* fft(stimArray,[],2), [], 2),2);
else
    blurArray = stimArray;
end

% Filter the spatially-blurred stimulus in time
if p.usePrefilter
    blurArray = filter(f.ex, 1, blurArray, [], 1);
end

lpArray = filter(f.lp, 1, blurArray, [], 1);
hpArray = filter(f.hp, 1, blurArray, [], 1);

% Shift the stimulus to get each of the photoreceptor inputs
% Note that signs of circshifts are reversed relative to index notation
prShift = floor(p.photoreceptorSpacing / p.dx);
resp1 = circshift(lpArray, +prShift, 2);
resp2 = hpArray;
resp3 = circshift(lpArray, -prShift, 2);

end

