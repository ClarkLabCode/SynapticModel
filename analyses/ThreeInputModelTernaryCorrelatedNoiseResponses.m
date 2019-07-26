function ThreeInputModelTernaryCorrelatedNoiseResponses(config, varargin)

%% Either load data from cache or run simulation

% Search for existing cache files
dataPath = dir(fullfile(config.rootDataPath, 'ternaryCorrelatedNoiseData_*.mat'));

% Set list of variable names to save or load
saveVarList = {'params','filters','noiseParam','tau','dtVec','meanResp','meanNumResp','meanDenResp'};

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
    
    % Delay times
    dtVec = [-12, -6:1:6, 12]';
    tau = dtVec / noiseParam.updateRate * 1000;
    
    % Allocate containers
    numT = length(dtVec);
    meanResp = nan(numT,3,params.numRep);
    meanNumResp = nan(numT,3,params.numRep);
    meanDenResp = nan(numT,3,params.numRep);
    
    % Iterate over realizations
    parfor indR = 1:params.numRep
        tic;
        
        % Iterate over correlation intervals
        for indT = 1:numT
            
            % Make the stimulus
            [ stimArray ] = TernaryCorrelatedNoise(params, noiseParam, dtVec(indT));
            
            % Compute the response
            [ meanResp(indT,:,indR), ~, ~, meanNumResp(indT,:,indR), meanDenResp(indT,:,indR) ] = ComputeThreeInputModelResponse(stimArray, params, filters);
            
        end
        
        % Print a status update
        fprintf('Realization %d of %d: %f seconds.\n', indR, params.numRep, toc);
    end
    
    % Save the data
    tic;
    savePath = fullfile(config.rootDataPath, sprintf('ternaryCorrelatedNoiseData_%s.mat', datestr(datetime('now'), 'yyyymmdd_HHMMSS')));
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

%% Define quantities for plotting

% Stimulus names
legendStr = {'positive correlations','negative correlations', 'uncorrelated'};

% Color order
corder = [flipud(lines(2)); 0.5,0.5,0.5];

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

%% Plot the raw responses

mu = nanmean(meanResp,3);
ci = bootci(params.nboot, {@nanmean, permute(meanResp,[3,1,2])});
cl = squeeze(ci(1,:,:,:));
cu = squeeze(ci(2,:,:,:));

MakeFigure;
PlotAsymmetricErrorPatch(tau, mu, cl, cu, corder);
legend(legendStr);
xlabel('interval (ms)');
ylabel('response (arb. units)');
% xlim([0 360]);
% xticks(0:45:360);
ylim([0 ceil(max(mu(:)))]);
ConfAxis(16);
axis('square');

%% Separate PD and ND
tauPD = tau(tau>=0);
muPD = mu(tau>=0,:);
muND = flipud(mu(tau<=0,:));
clPD = cl(tau>=0,:);
clND = flipud(cl(tau<=0,:));
cuPD = cu(tau>=0,:);
cuND = flipud(cu(tau<=0,:));

%% Plot positive correlations only

MakeFigure;
hold on;
PlotAsymmetricErrorPatch(tauPD, muPD(:,1), clPD(:,1),cuPD(:,1), corder(1,:), '-');
PlotAsymmetricErrorPatch(tauPD, muND(:,1), clND(:,1),cuND(:,1), corder(1,:), '--');
PlotAsymmetricErrorPatch(tauPD, muPD(:,3), clPD(:,3),cuPD(:,3), corder(3,:), '-');
legend({'PD','ND','UN'});
xlabel('interval (ms)');
ylabel('response (arb. units)');
ylim([0 ceil(max(mu(:)))]);
% xticks(0:45:360);
title('positive correlations');
ConfAxis(16);
axis('square');

%% Plot negative correlations only

MakeFigure;
hold on;
PlotAsymmetricErrorPatch(tauPD, muPD(:,2), clPD(:,2),cuPD(:,2), corder(2,:), '-');
PlotAsymmetricErrorPatch(tauPD, muND(:,2), clND(:,2),cuND(:,2), corder(2,:), '--');
PlotAsymmetricErrorPatch(tauPD, muPD(:,3), clPD(:,3),cuPD(:,3), corder(3,:), '-');
legend({'PD','ND','UN'});
xlabel('interval (ms)');
ylabel('response (arb. units)');
ylim([0 ceil(max(mu(:)))]);
% xticks(0:45:360);
title('negative correlations');
ConfAxis(16);
axis('square');

%% Plot the responses relative to uncorrelated noise

meanUn = mean(meanResp(:,3,:),3);
relResp = (meanResp - meanUn) ./ meanUn;
mu = nanmean(relResp,3);
ci = bootci(params.nboot, {@nanmean, permute(relResp,[3,1,2])});
cl = squeeze(ci(1,:,:,:));
cu = squeeze(ci(2,:,:,:));

MakeFigure;
PlotAsymmetricErrorPatch(tau, mu, cl,cu, corder);
legend(legendStr);
xlabel('interval (ms)');
ylabel('response relative to uncorrelated');
a = ceil(max(abs(relResp(:))));
ylim([-a,a]);
ConfAxis(16);
axis('square');

%% Plot the DSI

idxPD = (tau>=0);
idxND = (tau<=0);

dsi = (meanResp(idxPD,:,:) - meanResp(idxND,:,:)) ./ (meanResp(idxPD,:,:) + meanResp(idxND,:,:));
mu = nanmean(dsi,3);
ci = bootci(params.nboot, {@nanmean, permute(dsi,[3,1,2])});
cl = squeeze(ci(1,:,:,:));
cu = squeeze(ci(2,:,:,:));

MakeFigure;
PlotAsymmetricErrorPatch(tau(idxPD), mu, cl,cu, corder);
legend(legendStr);
xlabel('interval (ms)');
ylabel('DSI');
ylim([-1,1]);
ConfAxis(16);
axis('square');

%% Plot all DSIs together

dsiNum = (meanNumResp(idxPD,:,:) - meanNumResp(idxND,:,:)) ./ (meanNumResp(idxPD,:,:) + meanNumResp(idxND,:,:));
dsiDen = (meanDenResp(idxPD,:,:) - meanDenResp(idxND,:,:)) ./ (meanDenResp(idxPD,:,:) + meanDenResp(idxND,:,:));

dsiAll = cat(2, dsi(:,1:2,:), dsiNum(:,1:2,:), dsiDen(:,1:2,:));
mu = nanmean(dsiAll,3);
ci = bootci(params.nboot, {@nanmean, permute(dsiAll,[3,1,2])});
cl = squeeze(ci(1,:,:,:));
cu = squeeze(ci(2,:,:,:));

MakeFigure;
PlotAsymmetricErrorPatch(tau(idxPD), mu, cl,cu, corder2(1:6,:));
legend({'full model +','full model -','numerator +','numerator -','denominator +','denominator -'});
xlabel('interval (ms)');
ylabel('DSI');
ylim([-1,1]);
ConfAxis(16);
axis('square');

end