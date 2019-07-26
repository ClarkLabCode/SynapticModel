function ThreeInputModelSinusoidResponses(config, varargin)

%% Either load data from cache or run simulation

% Search for existing cache files
dataPath = dir(fullfile(config.rootDataPath, 'sinusoidBarPlotData_*.mat'));

% Set list of variable names to save or load
saveVarList = {'params','filters','c','tf','sf','stimArrayM','meanResp','voltageResp','calciumResp', 'meanNumResp','meanDenResp','meanNumRespLN'};

% Check if data regeneration flag is set or if no cache exists
if config.regenerateData || isempty(dataPath)
    fprintf('Simulation run flag set or no cache found. Running simulation...\n');
    tic;
    
    % Set overall parameters
    [ params ] = SetModelParameters('tInt', 0, varargin{:});
    
    % Make the filters
    [ filters ] = MakeModelFilters(params);
    
    % Contrast
    c = 1/2;
    
    % Temporal and spatial frequency vectors
    tf = 1;
    sf = 1/45;
    
    tic;
    
    % Compute monocomponent grating responses
    [ stimArrayM ] = MonocomponentSinusoidalGratings(params, c, tf, sf);
    [ meanRespM, voltageRespM, calciumRespM, meanNumRespM, meanDenRespM, meanNumRespLNM ] = ComputeThreeInputModelResponse(stimArrayM, params, filters);
    
    % Compute composite grating responses
    [ stimArrayC ] = CompositeSinusoidalGratings(params, c, tf, sf);
    [ meanRespC, voltageRespC, calciumRespC, meanNumRespC, meanDenRespC, meanNumRespLNC ] = ComputeThreeInputModelResponse(stimArrayC, params, filters);
    
    % Combine results together
    meanResp = cat(1, meanRespM, meanRespC);
    meanNumResp = cat(1, meanNumRespM, meanNumRespC);
    meanDenResp = cat(1, meanDenRespM, meanDenRespC);
    meanNumRespLN = cat(1, meanNumRespLNM, meanNumRespLNC);
    voltageResp = cat(3, voltageRespM, mean(voltageRespC,4));
    calciumResp = cat(3, calciumRespM, mean(calciumRespC,4));
    
    % Print a status update
    fprintf('Completed simulation in %f seconds\n', toc);
    
    % Save the data
    tic;
    savePath = fullfile(config.rootDataPath, sprintf('sinusoidBarPlotData_%s.mat', datestr(datetime('now'), 'yyyymmdd_HHMMSS')));
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

%% Set options for plotting

% Set colormap
load('utils/blueRedColorMap.mat', 'cmpRed','cmpBlueRed');

% Legend
legendStr = {'PD','ND','PD+ND','PD+OD'};

% Color order (from Badwan)
corder = [21 114 186; 216 85 39; 129 47 140; 60, 181, 74]/256;

%% Make bar plot of mean responses

MakeFigure;
bar(1:4, meanResp / meanResp(1), 'FaceColor','flat','CData', corder);
axis('square');
xticks(1:4);
xticklabels(legendStr);
ylabel('normalized response');
ConfAxis(16);

%% Make bar plot of mean numerator LNLN responses

MakeFigure;
bar(1:4, meanNumResp / meanNumResp(1), 'FaceColor','flat','CData', corder);
axis('square');
xticks(1:4);
xticklabels(legendStr);
ylabel('normalized response');
title('numerator LNLN');
ConfAxis(16);

%% Make bar plot of mean denominator LNLN responses

MakeFigure;
bar(1:4, meanDenResp / meanDenResp(1), 'FaceColor','flat','CData', corder);
axis('square');
xticks(1:4);
xticklabels(legendStr);
ylabel('normalized response');
title('denominator LNLN');
ConfAxis(16);

%% Make bar plot of mean numerator LN responses

MakeFigure;
bar(1:4, meanNumRespLN / meanNumRespLN(1), 'FaceColor','flat','CData', corder);
axis('square');
xticks(1:4);
xticklabels(legendStr);
ylabel('normalized response');
title('numerator LN');
ConfAxis(16);

end

