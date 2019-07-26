function [ f ] = MakeModelFilters( p )

%% Truncate temporal vector to positive times

t = p.t;
t = t(t>=0);
f.t = t;

%% Define temporal filters
% Note that these filters have unit L2 norm

if p.usePrefilter
    % Exponential photoreceptor filter
    f.ex = sqrt(p.dt) .* sqrt(2) .* (p.tauPr1^(-1/2)) .* double(t>=0) .* exp(-t/p.tauPr1);
    
    % Second order exponential lowpass and its derivative
    f.lp = sqrt(p.dt) .* 2 .* (p.tauPr2^(-3/2)) .* double(t>=0) .* t .* exp(-t/p.tauPr2);
    f.hp = sqrt(p.dt) .* 2 .* (p.tauPr2^(-3/2)) .* double(t>=0) .* (p.tauPr2 - t) .* exp(-t/p.tauPr2);
else
    
    % Second order exponential lowpass and its derivative
    f.lp = sqrt(p.dt) .* 2 .* (p.tauLp^(-3/2)) .* double(t>=0) .* t .* exp(-t/p.tauLp);
    f.hp = sqrt(p.dt) .* 2 .* (p.tauHp^(-3/2)) .* double(t>=0) .* (p.tauHp - t) .* exp(-t/p.tauHp);
end

%% Define spatial filters

% Convert FWHM to STD
sBlur = p.fwhmBlur / (2*sqrt(2*log(2)));

% Compute spatial filter
% Note that this filter has unit L1 norm
f.spatialFilter = exp(-(p.x - p.xExtent/2).^2/(2*sBlur^2)) / sqrt(2 * pi * sBlur^2) * p.dx;

end