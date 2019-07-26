function [ config ] = SetConfiguration(varargin)
%Function to set file paths

%% Set defaults

% Root data path to store cache files
config.rootDataPath = 'EDIT';

% Path to natural scene source data 
config.sceneSourcePath = 'EDIT';

% Select whether to regenerate data by default
config.regenerateData = false;

% Select whether to re-convert natural scenes to contrast by default
config.regenerateSceneData = false;

%% Extract user inputs

% Check for input arguments
if nargin > 0

    % Validate the number of input arguments
    if rem(nargin,2) ~= 0
        error('Arguments must be in name-value pairs!');
    end

    % Extract name-value pairs and store in parameter structure
    for ind = 1:2:nargin
        config.(varargin{ind}) = varargin{ind+1};
    end
end

end
