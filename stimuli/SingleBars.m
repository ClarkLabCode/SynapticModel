function [ stimArray ] = SingleBars(params, barParam)
%SINGLEBARS: A function to make single-bar stimuli. 

numStim = 2;

mlum = barParam.mlum;
c = barParam.c;
delay = barParam.delay;
barWidth = barParam.barWidth;
barOffset = barParam.barOffset;
barPeriod = barParam.barPeriod;

t = params.t;
x = params.x;
barOneMask = params.mask;

%% Make the individual bars

barOne   = +1 * c * barOneMask * barfun(x,barPeriod,barWidth);
barTwo   = -1 * c * barOneMask * barfun(x,barPeriod,barWidth);

%% Make the stimulus array

stimArray = cat(3, mlum + barOne, mlum + barTwo);

end

function [ bar ] = barfun(x, barPeriod, barWidth)
bar = double(mod(x, barPeriod) < barWidth);
end