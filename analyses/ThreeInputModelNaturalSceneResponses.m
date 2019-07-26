function ThreeInputModelNaturalSceneResponses(config, varargin)

%% Either load data from cache or run simulation

% Search for existing cache files
dataPath = dir(fullfile(config.rootDataPath, 'naturalSceneResponses_*.mat'));

% Set options for bootstrapping
opts = {'type','cper','options',statset('UseParallel',true)};

% Set list of variable names to save or load
saveVarList = {'params','filters','sceneDataPath','numScenes','v','rowSampleVec','meanHrcResp','meanRectHrcResp','meanResp','coactPerScene','meanRespAll','opts'};

% Check if data regeneration flag is set or if no cache exists
if config.regenerateData || config.regenerateSceneData || isempty(dataPath)
    
    % Set parameters
    [ params ] = SetModelParameters('tInt', 0,'useSpatialFilter', false, varargin{:});
    [ filters ] = MakeModelFilters(params);
    
    % Search for natural scene data file
    sceneDataPath = dir(fullfile(config.rootDataPath, 'naturalSceneData_*.mat'));
    
    % If no scene data cache exists or if the data regeneration flag is
    % set, regenerate natural scene data cache
    if config.regenerateSceneData || isempty(sceneDataPath)
        ConvertNaturalScenesToContrast(config, params);
        sceneDataPath = dir(fullfile(config.rootDataPath, 'naturalSceneData_*.mat'));
    end
    
    % Get the path to the most recent scene data cache
    sceneDataPath = fullfile(config.rootDataPath, sceneDataPath(end).name);
    
    % Load the natural scene data
    tic;
    load(sceneDataPath, 'contrastScenes');
    fprintf('Loaded natural scene data from %s in %f seconds\n', sceneDataPath, toc);
    
    % Get the size of the scene database
    [numY,numX,numScenes] = size(contrastScenes);
    
    % Get the spatial resolution of the scenes in degrees per pixel
    sceneResX = 360 / numX;
    
    % Define scene x and y-vectors
    sceneVecX = (0:numX-1)' * sceneResX;
    
    % Convert photoreceptor spacing to pixels
    prSpacingPx = params.photoreceptorSpacing / sceneResX;
    
    % Sample velocities from a Gaussian distribution with enforced symmetry
    v = abs(randn(params.numRep/2,1)) * params.velStd;
    v = [-v; v];
    
    % Draw scene samples IID from a uniform distribution with enforced symmetry
    sceneSampleVec = repmat(randi([1,numScenes], params.numRep/2,1),2,1);
    
    % Sample rows IID from a uniform distribution with enforced symmetry
    rowSampleVec = repmat(randi([1, numY], params.numRep/2, 1),2,1);
    
    % Sample columns IID from a uniform distribution with enforced symmetry
    colSampleVec = repmat(randi([1, numX], params.numRep/2, 1),2,1);
    
    % Allocate containers
    numV = length(v);
    meanHrcResp = nan(numV,1);
    meanRectHrcResp = nan(numV,1);
    meanResp = nan(numV,4);
    coactPerScene = nan(numV,4,4);
    
    % Iterate over velocities
    parfor indV = 1:numV
        tic;
        
        % Compute the offsets at each timestep (in pixels)
        dx = -params.t * v(indV) / sceneResX + colSampleVec(indV);
        
        % Compute the indices for each input with circular boundary
        % conditions
        idx1 = 1 + mod(round(dx - prSpacingPx), numX);
        idx2 = 1 + mod(round(dx), numX);
        idx3 = 1 + mod(round(dx + prSpacingPx), numX);
        
        % Extract the desired row
        sampleScene = squeeze(contrastScenes(rowSampleVec(indV),:,sceneSampleVec(indV)))';
        
        % Filter inputs in time
        resp1 = filter(filters.lp, 1, sampleScene(idx1), [], 1);
        resp2 = filter(filters.hp, 1, sampleScene(idx2), [], 1);
        resp3 = filter(filters.lp, 1, sampleScene(idx3), [], 1);
        resp4 = filter(filters.hp, 1, sampleScene(idx1), [], 1);
        resp5 = filter(filters.lp, 1, sampleScene(idx2), [], 1);
        
        % Compute HRC response
        hrcResp = resp1 .* resp2 - resp4 .* resp5;
        relu = @(f) (f .* (f>0));
        [ meanHrcResp(indV) ] = mean(hrcResp(params.averagingMask));
        [ meanRectHrcResp(indV) ] = mean(relu(hrcResp(params.averagingMask)));
        
        % Compute model response
        [ calciumResp ] = ComputeThreeInputModelSplitChannelResponsesFromFilteredInputs(resp1, resp2, resp3, params);
        calciumResp = squeeze(calciumResp(params.averagingMask,:,:,:,:,:));
        meanResp(indV,:) = mean(calciumResp,1);
        
        % Compute coactivations
        coactPerScene(indV,:,:) = calciumResp' * calciumResp;
        
        % Print a status update
        if ~mod(indV, 1e4)
            fprintf('Velocity %d of %d: %f seconds.\n', indV, numV, toc);
        end
    end
    
    % Concatenate all responses into one array
    meanRespAll = cat(2, meanHrcResp, meanRectHrcResp, meanResp);
    
    % Save the data
    tic;
    savePath = fullfile(config.rootDataPath, sprintf('naturalSceneResponses_%s.mat', datestr(datetime('now'), 'yyyymmdd_HHMMSS')));
    save(savePath, '-v7.3', saveVarList{:});
    fprintf('Saved data to %s in %f seconds\n', savePath, toc);
    
else
    
    % Get the path to the most recent file
    dataPath = fullfile(config.rootDataPath, dataPath(end).name);
    
    % Load the data
    tic;
    load(dataPath, saveVarList{:});
    fprintf('Loaded data from %s in %f seconds\n', dataPath, toc);
    
end

%% Define quantities for plotting

% Color order for bar plots
corder = repmat(lines(1), 6, 1);

% Other color order
corder2 = lines(6);

cellLabelStr = {'T4p','T4r','T5p','T5r'};

load('utils/blueRedColorMap.mat','cmpRed','cmpBlueRed');

%% Plot mean response as a function of velocity

binEdges = (-300:10:300)';
binCenters = binEdges(1:end-1) + diff(binEdges)/2;
vDisc = discretize(v, binEdges);
notNanIdx = ~isnan(vDisc);

mu = nan(length(binCenters), size(meanRespAll,2));
ci = nan(length(binCenters), size(meanRespAll,2), 2);
for ind = 1:size(meanRespAll,2)
    mu(:,ind) = accumarray(vDisc(notNanIdx), meanRespAll(notNanIdx,ind), [length(binCenters),1], @mean);
    ci(:,ind,:) = permute(bootci(params.nboot, {@(x,y) accumarray(x,y, [length(binCenters),1], @mean), vDisc(notNanIdx), meanRespAll(notNanIdx,ind)}, opts{:}), [2,3,1]);
end

MakeFigure;
histogram(v, binEdges, 'normalization','pdf','edgecolor','none');
xlabel('velocity (\circ/s)');
ylabel('pdf (s/\circ)');
axis('square');
ConfAxis(16);

MakeFigure;
PlotAsymmetricErrorPatch(binCenters, mu(:,1:2), ci(:,1:2,1), ci(:,1:2,2), corder2(1:2,:));
xlabel('velocity (\circ/s)');
ylabel('response (arb. units)');
legend({'HRC', 'HRCr'});
axis('square');
ConfAxis(16);

MakeFigure;
PlotAsymmetricErrorPatch(binCenters, mu(:,3:end), ci(:,3:end,1), ci(:,3:end,2), corder2(3:end,:));
xlabel('velocity (\circ/s)');
ylabel('response (arb. units)');
legend(cellLabelStr);
axis('square');
ConfAxis(16);

%% Show the coactivations

% Average the coactivations over scenes and velocities
rawCoact = squeeze(nanmean(coactPerScene,1));
MakeFigure;
imagesc(1:4, 1:4, rawCoact ./ diag(rawCoact));
axis('xy','square');
xticks(1:4);
xticklabels(cellLabelStr);
yticks(1:4);
yticklabels(cellLabelStr);
colormap(cmpRed);
caxis([0 1]);
cbar = colorbar;
cbar.Ticks = 0:0.5:1;
ylabel(cbar, 'coactivation (arb. units)');
ConfAxis(16);

end