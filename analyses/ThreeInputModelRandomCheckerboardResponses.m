function ThreeInputModelRandomCheckerboardResponses(config, varargin)


%% Either load data from cache or run simulation

% Search for existing cache files
dataPath = dir(fullfile(config.rootDataPath, 'randomCheckerboardData_*.mat'));

% Set list of variable names to save or load
saveVarList = {'params','filters','noiseParam','v','meanResp','meanNumResp','meanDenResp'};

% Check if data regeneration flag is set or if no cache exists
if config.regenerateData || isempty(dataPath)
    fprintf('Simulation run flag set or no cache found. Running simulation...\n');
    tic;

    % Set overall parameters
    [ params ] = SetModelParameters('tInt', 0, varargin{:});
    [ filters ] = MakeModelFilters(params);

    % General noise parameters
    noiseParam.barWidth = 5;
    noiseParam.updateRate = 60;
    noiseParam.dx = 1;
    noiseParam.c = 1;
    noiseParam.mlum = 0;
    noiseParam.density = 0.4;

    % Velocities
    v = (100)';
    numV = length(v);

    % Allocate containers
    meanResp = nan(numV,6,params.numRep);
    meanNumResp = nan(numV,6,params.numRep);
    meanDenResp = nan(numV,6,params.numRep);

    % Iterate over realizations
    parfor indR = 1:params.numRep
        tic;

        % Iterate over velocities
        for indV = 1:numV

            % Make the stimulus
            [ stimArray ] = RandomCheckerboard(params, noiseParam, v(indV));

            % Compute the response
            [ meanResp(indV,:,indR), ~, ~, meanNumResp(indV,:,indR), meanDenResp(indV,:,indR) ] = ComputeThreeInputModelResponse(stimArray, params, filters);

        end

        % Print a status update
        fprintf('Realization %d of %d: %f seconds.\n', indR, params.numRep, toc);
    end

    % Save the data
    tic;
    savePath = fullfile(config.rootDataPath, sprintf('randomCheckerboardData_%s.mat', datestr(datetime('now'), 'yyyymmdd_HHMMSS')));
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

%% Set preferences for plotting

% Stimulus names
legendStr = {'PD','ND','PD+ND','PD+OD','ND+OD','OD'};

% Color order
corder = [lines(5);0.5,0.5,0.5];

corder2 = [
    0.1059    0.6196    0.4667
    0.8510    0.3725    0.0078
    0.4588    0.4392    0.7020
    0.9059    0.1608    0.5412
    0.6549    0.3961    0.3176
    0.4000    0.6510    0.1176
    0.9020    0.6706    0.0078
    0.6510    0.4627    0.1137
    0.4000    0.4000    0.4000
    ];

corder3 = [21 114 186; 216 85 39; 129 47 140; 60, 181, 74; 50+60, 50+181, 50+74]/256;

%% Plot the raw responses

% Compute statistics
mu = nanmean(meanResp(:,1:5,:),3)';
ci = bootci(params.nboot, {@nanmean, permute(meanResp(:,1:5,:),[3,1,2])});
cl = squeeze(ci(1,:,:,:));
cu = squeeze(ci(2,:,:,:));

% Make the bar plot
MakeFigure;
ErrorBarChart((1:5)', mu, cl, cu, corder3);
xticks(1:5);
xticklabels(legendStr(1:5));
ylabel('response (arb. units)');
title(sprintf('%d \\circ/s', v));
axis('square');
ConfAxis(16);

%% Plot the raw numerator and denominator LNLN responses

mu = nanmean(meanNumResp(:,1:5,:),3)';
ci = bootci(params.nboot, {@nanmean, permute(meanNumResp(:,1:5,:),[3,1,2])});
cl = squeeze(ci(1,:,:,:));
cu = squeeze(ci(2,:,:,:));

MakeFigure;
ErrorBarChart((1:5)', mu, cl, cu, corder3);
xticks(1:5);
xticklabels(legendStr(1:5));
ylabel('numerator LNLN response (arb. units)');
title(sprintf('%d \\circ/s', v));
axis('square');
ConfAxis(16);
title('Numerator LNLN')

mu = nanmean(meanDenResp(:,1:5,:),3)';
ci = bootci(params.nboot, {@nanmean, permute(meanDenResp(:,1:5,:),[3,1,2])});
cl = squeeze(ci(1,:,:,:));
cu = squeeze(ci(2,:,:,:));

MakeFigure;
ErrorBarChart((1:5)', mu, cl, cu, corder3);
xticks(1:5);
xticklabels(legendStr(1:5));
ylabel('denominator LNLN response (arb. units)');
title(sprintf('%d \\circ/s', v));
axis('square');
ConfAxis(16);
title('Denominator LNLN')

end
