function ThreeInputModelOrientedGratingResponses(config, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Either load data from cache or run simulation

% Search for existing cache files
dataPath = dir(fullfile(config.rootDataPath, 'orientedGratingData_*.mat'));

% Set list of variable names to save or load
saveVarList = {'params','filters','c', 'tf','sf','theta', 'meanResp'};

% Check if data regeneration flag is set or if no cache exists
if config.regenerateData || isempty(dataPath)
    fprintf('Simulation run flag set or no cache found. Running simulation...\n');
    tic;
    
    % Set overall parameters
    [ params ] = SetModelParameters('tInt', 0, 'dx', 0.5, varargin{:});
    
    % Make the filters
    [ filters ] = MakeModelFilters(params);
    
    % Contrast
    c = 1/2;
    
    % Temporal frequency
    tf = [0;1];
    
    % Spatial frequency
    sf = 1/45;
    
    % Orientation
    theta = (0:5:360)';
    
    % Compute responses at each TF and orientation
    numT = length(tf);
    numO = length(theta);
    meanResp = nan(numO,numT);
    for indT = 1:numT
        parfor indO = 1:numO
            tic;
            
            [ stimArray ] = OrientedSinusoidGratings(params, filters, c, tf(indT), sf, deg2rad(theta(indO)));
            
            [ meanResp(indO,indT) ] = ComputeThreeInputModelResponse(stimArray, params, filters);
            
            fprintf('TF %d of %d, orientation %d of %d: %f s.\n', indT, numT, indO, numO,toc);
        end
    end
    
    % Save the data
    tic;
    savePath = fullfile(config.rootDataPath, sprintf('orientedGratingData_%s.mat', datestr(datetime('now'), 'yyyymmdd_HHMMSS')));
    save(savePath, '-v7.3', saveVarList{:});
    fprintf('Saved data to %s in %f seconds\n', savePath, toc);
    
else
    tic;
    
    % Get the path to the most recent file
    dataPath = fullfile(config.rootDataPath, dataPath(end).name);
    
    % Load the data
    load(dataPath, saveVarList{:});
    
    % Print a status update
    fprintf('Loaded data from %s in %f seconds\n', dataPath, toc);
    
end


%% Plot the responses

% Symmetrize static responses to mod 180 as in Fisher
tMod180 = unique(mod(theta, 180));
meanRespSym = meanResp;
for ind = 1:length(tMod180)
    idx = mod(theta, 180) == tMod180(ind);
   meanRespSym(idx,1) = mean(meanResp(idx,1));
end


MakeFigure;
p = polaraxes();
polarplot(deg2rad(theta), meanRespSym ./ meanRespSym(1,:), '-o', 'linewidth', 2);
legend({'static', '1 Hz'});
p.FontSize = 16;
p.RAxis.Label.String = 'normalized response';
p.RLim = [0 1.25];
p.ThetaAxis.Label.String = 'orientation (\circ)';

end

