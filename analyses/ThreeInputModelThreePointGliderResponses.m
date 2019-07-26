function ThreeInputModelThreePointGliderResponses(config, varargin)

%% Either load data from cache or run simulation

% Search for existing cache files
dataPath = dir(fullfile(config.rootDataPath, 'threePointGliderData_*.mat'));

% Set list of variable names to save or load
% saveVarList = {'params','filters','noiseParam','updateRate','meanResp','meanNumResp','meanDenResp'};
saveVarList = {'params','filters','noiseParam','updateRate','meanResp'};

% Check if data regeneration flag is set or if no cache exists
if config.regenerateData || isempty(dataPath)
    fprintf('Simulation run flag set or no cache found. Running simulation...\n');
    tic;
    
    % Set overall parameters
    [ params ] = SetModelParameters('tInt', 0, varargin{:});
    [ filters ] = MakeModelFilters(params);
    
    % General glider parameters
    noiseParam.barWidth = 5;
    noiseParam.dx = 1;
    noiseParam.c = 1/2;
    noiseParam.mlum = 0;
    
    % Update rates
    updateRate = [(0.5:0.5:10)'; (12:2:60)'];
    
    % Get the number of update rates
    numT = length(updateRate);
    
    % Allocate a container
    meanResp = nan(numT,9,params.numRep);
    %     meanNumResp = nan(numT,9,params.numRep);
    %     meanDenResp = nan(numT,9,params.numRep);
    
    % Iterate over realizations
    parfor indR = 1:params.numRep
        
        % Start a timer for this realization
        tic;
        
        % Iterate over update rates
        for indT = 1:numT
            
            % Make the stimulus
            [ stimArray ] = ThreePointGlider(params, noiseParam, updateRate(indT));
            
            % Compute the response
            %             [ meanResp(indT,:,indR), ~, ~, meanNumResp(indT,:,indR), meanDenResp(indT,:,indR) ] = ComputeThreeInputModelResponse(stimArray, params, filters);
            [ meanResp(indT,:,indR) ] = ComputeThreeInputModelResponse(stimArray, params, filters);
            
        end
        
        % Print a status update
        fprintf('Realization %d of %d: %f seconds.\n', indR, params.numRep, toc);
    end
    
    
    % Save the data
    tic;
    savePath = fullfile(config.rootDataPath, sprintf('threePointGliderData_%s.mat', datestr(datetime('now'), 'yyyymmdd_HHMMSS')));
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
legendStr = {'PDcon+','PDcon-','NDcon+','NDcon-','PDdiv+','PDdiv-','NDdiv+','NDdiv-','UN'};

% Color order
corder = [
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

corder2 = lines(4);

%% Plot all raw responses

mu = nanmean(meanResp,3);
ci = bootci(params.nboot, {@nanmean, permute(meanResp,[3,1,2])});
cl = squeeze(ci(1,:,:,:));
cu = squeeze(ci(2,:,:,:));

if length(updateRate) < 6
    for ind = 1:length(updateRate)
        MakeFigure;
        ErrorBarChart(1:size(mu,2), mu(ind,:), cl(ind,:),cu(ind,:), corder);
        xticks(1:size(mu,2));
        xticklabels(legendStr);
        ylabel('response (arb. units)');
        title(sprintf('%d Hz update rate',updateRate(ind)));
        ConfAxis(16);
        axis('square');
    end
else
    MakeFigure;
    PlotAsymmetricErrorPatch(updateRate, mu, cl,cu, corder);
    xlabel('update rate (Hz)');
    legend(legendStr);
    ylabel('response (arb. units)');
    ConfAxis(16);
    axis('square');
end

%% Plot all responses relative to uncorrelated

meanUn = mean(meanResp(:,end,:),3);
relResp = (meanResp(:, 1:end-1, :) - meanUn) ./ meanUn;
mu = nanmean(relResp,3);
ci = bootci(params.nboot, {@nanmean, permute(relResp,[3,1,2])});
cl = squeeze(ci(1,:,:,:));
cu = squeeze(ci(2,:,:,:));

if length(updateRate) < 6
    for ind = 1:length(updateRate)
        MakeFigure;
        ErrorBarChart(1:size(mu,2), mu(ind,:), cl(ind,:),cu(ind,:), corder(1:end-1,:));
        xticks(1:size(mu,2));
        xticklabels(legendStr);
        ylabel('response relative to uncorrelated');
        title(sprintf('%d Hz update rate',updateRate(ind)));
        ConfAxis(16);
        axis('square');
    end
else
    MakeFigure;
    PlotAsymmetricErrorPatch(updateRate, mu, cl,cu, corder(1:end-1,:));
    hold on;
    plot([0 max(updateRate)], [0 0], '--k', 'linewidth', 2);
    xlabel('update rate (Hz)');
    legend(legendStr);
    ylabel('response relative to uncorrelated');
    ConfAxis(16);
    axis('square');
end

%% Plot differences of PD and ND

respDiff = meanResp(:,contains(legendStr, 'PD'),:) - meanResp(:,contains(legendStr, 'ND'),:);
mu = nanmean(respDiff,3);
ci = bootci(params.nboot, {@nanmean, permute(respDiff,[3,1,2])});
cl = squeeze(ci(1,:,:,:));
cu = squeeze(ci(2,:,:,:));

if length(updateRate) < 6
    for ind = 1:length(updateRate)
        MakeFigure;
        ErrorBarChart(1:size(mu,2), mu(ind,:), cl(ind,:),cu(ind,:), corder2);
        xticks(1:size(mu,2));
        xticklabels(erase(erase(legendStr([1,2,5,6]),'PD'), 'ND'));
        ylabel('net response (arb. units)');
        title(sprintf('%d Hz update rate',updateRate(ind)));
        ConfAxis(16);
        axis('square');
    end
else
    MakeFigure;
    PlotAsymmetricErrorPatch(updateRate, mu, cl,cu, corder2);
    hold on;
    plot([0 max(updateRate)], [0 0], '--k', 'linewidth', 2);
    xlabel('update rate (Hz)');
    legend(erase(erase(legendStr([1,2,5,6]),'PD'), 'ND'));
    ylabel('net response (arb. units)');
    ConfAxis(16);
    axis('square');
end


%% Plot bar graphs of net responses

updateRateList = [5;10;24;40;60];

for ind = 1:length(updateRateList)
    
    % Extract the desired responses
    respSel = respDiff(updateRate == updateRateList(ind),:,:);
    
    % Compute statistics
    mu = nanmean(respSel,3)';
    ci = bootci(params.nboot, {@nanmean, permute(respSel,[3,1,2])});
    cl = squeeze(ci(1,:,:,:));
    cu = squeeze(ci(2,:,:,:));
    
    % Make the bar plot
    MakeFigure;
    ErrorBarChart((1:4)', mu, cl, cu, corder2);
    xticks(1:4);
    xticklabels(erase(erase(legendStr([1,2,5,6]),'PD'), 'ND'));
    ylabel('net response (arb. units)');
    title(sprintf('%d Hz', updateRateList(ind)));
    axis('square');
    ConfAxis(16);
    
end

end