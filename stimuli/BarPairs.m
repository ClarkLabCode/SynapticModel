function [ stimArray ] = BarPairs(params, barParam)
% Make all bar direction + polarity combinations

numStim = 8;

stimArray = nan(length(params.t),length(params.x),numStim);

% ON-ON
stimArray(:,:,1:2) = MakeOneBarPair(params, barParam, +1, +1);

% OFF-OFF
stimArray(:,:,3:4) = MakeOneBarPair(params, barParam, -1, -1);

% ON-OFF
stimArray(:,:,5:6) = MakeOneBarPair(params, barParam, +1, -1);

% OFF-ON
stimArray(:,:,7:8) = MakeOneBarPair(params, barParam, -1, +1);

end


function [ stimArray ] = MakeOneBarPair(params, barParam, barOnePol, barTwoPol)
% Make one PD-ND bar pairing

mlum = barParam.mlum;
c = barParam.c;
delay = barParam.delay;
barWidth = barParam.barWidth;
barOffset = barParam.barOffset;
barPeriod = barParam.barPeriod;

t = params.t;
x = params.x;
barOneMask = params.mask;
barTwoMask = params.mask & (t >= delay);

%% Make the individual bars 
% (locations are define relative to the lagging bar as in Salazar-Gatizmas et al. 2018)

barTwo   = c * barTwoPol * barTwoMask * barfun(x,barPeriod,barWidth);
barOnePD = c * barOnePol * barOneMask * barfun(x+barOffset,barPeriod,barWidth);
barOneND = c * barOnePol * barOneMask * barfun(x-barOffset,barPeriod,barWidth);

%% Make the stimulus array

stimArray = cat(3, mlum + barTwo + barOnePD, mlum + barTwo + barOneND);

end

function [ bar ] = barfun(x, barPeriod, barWidth)
bar = double(mod(x, barPeriod) < barWidth);
end