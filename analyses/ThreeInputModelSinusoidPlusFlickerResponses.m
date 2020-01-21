function ThreeInputModelSinusoidPlusFlickerResponses(config, tf, sf, varargin)

%% Either load data from cache or run simulation

% Check arguments
if nargin < 2, tf = 1; end
if nargin < 3, sf = 1/45; end

% Search for existing cache files
dataPath = dir(fullfile(config.rootDataPath, 'sinusoidPlusFlickerData_*.mat'));

% Set list of variable names to save or load
saveVarList = {'params','filters','noiseParam','tf','sf','cNoise','cSin','meanResp','meanHrcResp'};

% Check if data regeneration flag is set or if no cache exists
if config.regenerateData || isempty(dataPath)
    fprintf('Simulation run flag set or no cache found. Running simulation...\n');
    tic;
    
    % Set overall parameters
    [ params ] = SetModelParameters('tInt', 0, 'dx', 0.5, 'tOn', 3, varargin{:});
    [ filters ] = MakeModelFilters(params);
    
    % General noise parameters
    noiseParam.updateRate = 30;
    noiseParam.dx = 1;
    noiseParam.c = 1;
    noiseParam.mlum = 0;
    
    % Contrasts
    cSin = 0.4;
    cNoise = [0; 0.1; 0.2; 0.4; 0.6];
    
    % Make sinusoid stimulus array
    [ sinusoidStimArray ] = MonocomponentSinusoidalGratings(params, cSin, tf, sf);
    
    % Iterate over realizations in parallel
    numC = length(cNoise);
    meanResp = nan(params.numRep, 2, numC);
    meanHrcResp = nan(params.numRep, 2, numC);
    
    for indC = 1:numC
        parfor indR = 1:params.numRep
            tic;
            % Make noise stimulus array
            [ noiseStimArray ] = FullFieldBinaryFlicker(params, noiseParam, cNoise(indC));
            
            % Sum stimuli
            stimArray = sinusoidStimArray + noiseStimArray
            
            % Compute response
            [ meanResp(indR,:,indC) ] = ComputeThreeInputModelResponse(stimArray, params, filters);
            
            % Compute HRC response
            [ meanHrcResp(indR,:,indC) ] = ComputeRectifiedHRCResponse(stimArray, params, filters);
            
            if ~mod(indR,100)
                fprintf('Noise level %d of %d, realization %d of %d: %f s.\n', indC, numC, indR, params.numRep, toc);
            end
        end
    end
    
    % Save the data
    tic;
    savePath = fullfile(config.rootDataPath, sprintf('sinusoidPlusFlickerData_%s.mat', datestr(datetime('now'), 'yyyymmdd_HHMMSS')));
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

corder = lines(4);
% legendStr = {'HRC PD','HRC ND','Synaptic PD','Synaptic ND'};
legendStr = {'[HRC]_+ PD','[HRC]_+ ND','Synaptic PD','Synaptic ND'};

%% Plot responses at noise constrast 0.6 relative to noiseless responses

idx = (cNoise == 0.6);
relResp = cat(2, meanHrcResp(:,:,idx) ./ meanHrcResp(:,1,1), meanResp(:,:,idx) ./ meanResp(:,1,1));
mu = mean(relResp,1);
ci = bootci(params.nboot, {@mean, relResp});

MakeFigure;
ErrorBarChart(1:4, mu, ci(1,:), ci(2,:), corder);
xticks(1:4);
xlim([0 5]);
xticklabels(legendStr);
ylabel('response relative to noise-free PD (arb. units)');
axis('square');
% ylim([-1.2, 1.2]);
ConfAxis(16);
title('Contrast 0.6 flicker, contrast 0.4 sinusoid at 1 Hz, 45\circ');

%% Plot sweeps

muHrc = squeeze(mean(meanHrcResp./meanHrcResp(:,1,1),1))';
ciHrc = permute(squeeze(bootci(params.nboot, {@mean, meanHrcResp./meanHrcResp(:,1,1)})), [3,2,1]);

muSyn = squeeze(mean(meanResp./meanResp(:,1,1),1))';
ciSyn = permute(squeeze(bootci(params.nboot, {@mean, meanResp./meanResp(:,1,1)})), [3,2,1]);

mu = cat(2, muHrc, muSyn);
ci = cat(2, ciHrc, ciSyn);

MakeFigure;
PlotAsymmetricErrorPatch(cNoise, mu, ci(:,:,1), ci(:,:,2), corder);
xlabel('flicker contrast');
ylabel('normalized response');
legend(legendStr, 'location','southwest');
axis('square');
ConfAxis(16);


end