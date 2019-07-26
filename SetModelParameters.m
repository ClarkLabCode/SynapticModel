function [ p ] = SetModelParameters(varargin)

%% Set timing parameters

% Interleave duration (used to pad beginning and end of stimulus)
p.tInt = 0;

% Stimulus duration
p.tOn = 3;

% Duration to exclude onset transients
p.tAv = 1;

% Temporal resolution (s)
p.dt = 1/240;

%% Set spatial parameters

% Spatial extent (fills entire 360-degree space)
p.xExtent = 360;

p.xPad = 30;

% Spatial resolution (degrees)
p.dx = 0.5;

%% Set parameters for spatial filters

p.useSpatialFilter = true;

% Blur filter FWHM (degrees)
p.fwhmBlur = 5.7;

% Averaging filter FWHM (degrees)
p.fwhmAverage = 20 * (2*sqrt(2*log(2)));

% Set photoreceptor spacing (degrees)
p.photoreceptorSpacing = 5;
% p.photoreceptorSpacing = 5.1;

%% Set temporal filter parameters

% Set whether to use exponential or Butterworth filters
p.usePrefilter = false;

% 
p.tauPr1 = 0.01;
p.tauPr2 = 0.02;

% Filter time constants (s)
p.tauLp = 0.15;
p.tauHp = 0.15;

%% Set spatial phase shifts
% Note that these phase shifts are not used in the linearity analysis,
% where the phase shifts are fixed to be (0:1:7)/8*pi as in Wienecke et al.

p.useRandomShifts = false;

%% Set hardness of rectifiers

% For soft rectification:
% params.inputRectBeta = 2;
% params.outputRectBeta = 32;

% For hard rectification:
p.inputRectBeta = Inf;
p.outputRectBeta = Inf;

%% Set reversal potentials & conductances

% Input 1 ~ Mi9
% Input 2 ~ Mi1
% Input 3 ~ Mi4
p.V1 = - 30;
p.V2 = + 60;
p.V3 = - 30;
p.gleak = 1;

p.g1 = 0.3;
p.g2 = 0.1;
p.g3 = 0.3;

%% Simulation parameters

% Number of realizations
p.numRep = 1000;

% Number of bootstraps
p.nboot = 1000;

% Batch size
p.batchSize = 100;

% Kernel extraction method
p.olsKernel = false;

%% Natural scene simulation parameters

% Standard deviation of Gaussian distribution of velocities
p.velStd = 100;

%% Extract user inputs

% Check for input arguments
if nargin > 0
    
    % Validate the number of input arguments
    if rem(nargin,2) ~= 0
        error('Arguments must be in name-value pairs!');
    end
    
    % Extract name-value pairs and store in parameter structure
    for ind = 1:2:nargin
        p.(varargin{ind}) = varargin{ind+1};
    end
end

%% Make temporal + spatial vectors and masks for stimulus generation

% Total simulation time
p.tTot = 2*p.tInt + p.tOn;

% Compute spatial position vector
p.xTot = p.xExtent;
p.x = (0:p.dx:p.xExtent-p.dx);

% Compute temporal position vector
if p.tInt <= 0
    p.t = (0:p.dt:p.tOn-p.dt)';
else
    p.t = (-p.tInt:p.dt:p.tOn+p.tInt-p.dt)';
end

% Define stimulus presentation mask
if p.tInt <= 0
    p.mask = ones(size(p.t), 'logical');
else
    p.mask = (p.t>=0) & (p.t<p.tOn);
end

% Define averaging mask
if p.tAv <= 0
    p.averagingMask = p.mask;
else
    p.averagingMask = (p.t >= p.tAv) & (p.t<p.tOn);
end

end

