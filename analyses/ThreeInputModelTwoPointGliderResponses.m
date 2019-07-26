function ThreeInputModelTwoPointGliderResponses(config, varargin)

%% Either load data from cache or run simulation

% Search for existing cache files
dataPath = dir(fullfile(config.rootDataPath, 'twoPointGliderData_*.mat'));

% Set list of variable names to save or load
saveVarList = {'params','filters','noiseParam','updateRate','meanResp','meanNumResp','meanDenResp'};

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
    meanResp = nan(numT,5,params.numRep);
    meanNumResp = nan(numT,5,params.numRep);
    meanDenResp = nan(numT,5,params.numRep);
    
    % Iterate over realizations
    parfor indR = 1:params.numRep
        
        % Start a timer for this realization
        tic;
        
        % Iterate over update rates
        for indT = 1:numT
            
            % Make the stimulus
            [ stimArray ] = TwoPointGlider(params, noiseParam, updateRate(indT));
            
            % Compute the response
            [ meanResp(indT,:,indR), ~, ~, meanNumResp(indT,:,indR), meanDenResp(indT,:,indR) ] = ComputeThreeInputModelResponse(stimArray, params, filters);
            
        end
        
        % Print a status update
        fprintf('Realization %d of %d: %f seconds.\n', indR, params.numRep, toc);
    end
    
    % Save the data
    tic;
    savePath = fullfile(config.rootDataPath, sprintf('twoPointGliderData_%s.mat', datestr(datetime('now'), 'yyyymmdd_HHMMSS')));
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
legendStr = {'PD+','ND+','PD-','ND-','UN'};

% Color order
corder = [lines(4); 0.5, 0.5, 0.5];

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

% Color order (from Badwan)
% corder3 = [21 114 186; 216 85 39; 129 47 140; 60, 181, 74]/256;
corder3 = [128, 128, 128; 21 114 186; 216 85 39]/256;

%% Plot bar graphs of responses at the update rates presented in Badwan

updateRateList = [5;20;60];

for ind = 1:length(updateRateList)
    
    % Extract the desired responses
    respSel = meanResp(updateRate == updateRateList(ind), [5,1,2],:);
    
    % Compute statistics
    mu = nanmean(respSel,3)';
    ci = bootci(params.nboot, {@nanmean, permute(respSel,[3,1,2])});
    cl = squeeze(ci(1,:,:,:));
    cu = squeeze(ci(2,:,:,:));
    
    % Make the bar plot
    MakeFigure;
    ErrorBarChart((1:3)', mu, cl, cu, corder3);
    xticks(1:3);
    xticklabels(legendStr([5,1,2]));
    ylabel('response (arb. units)');
    title(sprintf('%d Hz', updateRateList(ind)));
    axis('square');
    ConfAxis(16);
    
end

%% Plot the raw responses

mu = nanmean(meanResp,3);
ci = bootci(params.nboot, {@nanmean, permute(meanResp,[3,1,2])});
cl = squeeze(ci(1,:,:,:));
cu = squeeze(ci(2,:,:,:));

MakeFigure;
PlotAsymmetricErrorPatch(updateRate, mu, cl,cu, corder);
legend(legendStr);
xlabel('update rate (Hz)');
ylabel('response (arb. units)');
ConfAxis(16);
axis('square');

%% Plot the responses relative to uncorrelated

meanUn = mean(meanResp(:,5,:),3);
relResp = (meanResp - meanUn) ./ meanUn;
mu = nanmean(relResp,3);
ci = bootci(params.nboot, {@nanmean, permute(relResp,[3,1,2])});
cl = squeeze(ci(1,:,:,:));
cu = squeeze(ci(2,:,:,:));

MakeFigure;
PlotAsymmetricErrorPatch(updateRate, mu, cl,cu, corder);
legend(legendStr);
xlabel('update rate (Hz)');
ylabel('response relative to uncorrelated');
ConfAxis(16);
axis('square');
a = ceil(max(abs(relResp(:))));
ylim([-a,a]);

%% Plot the DSI

dsi = (meanResp(:,[1,3],:)-meanResp(:,[2,4],:)) ./ (meanResp(:,[1,3],:)+meanResp(:,[2,4],:));
mu = nanmean(dsi,3);
ci = bootci(params.nboot, {@nanmean, permute(dsi,[3,1,2])});
cl = squeeze(ci(1,:,:,:));
cu = squeeze(ci(2,:,:,:));

MakeFigure;
PlotAsymmetricErrorPatch(updateRate, mu, cl,cu, corder);
legend({'+','-'});
xlabel('update rate (Hz)');
ylabel('DSI');
ConfAxis(16);
axis('square');
ylim([-1,1]);

%% Plot the raw numerator and denominator LNLN responses

mu = nanmean(meanNumResp,3);
ci = bootci(params.nboot, {@nanmean, permute(meanNumResp,[3,1,2])});
cl = squeeze(ci(1,:,:,:));
cu = squeeze(ci(2,:,:,:));

MakeFigure;
PlotAsymmetricErrorPatch(updateRate, mu, cl,cu, corder);
legend(legendStr);
xlabel('update rate (Hz)');
ylabel('numerator LNLN response (arb. units)');
ConfAxis(16);
title('Numerator LNLN');
axis('square');

mu = nanmean(meanDenResp,3);
ci = bootci(params.nboot, {@nanmean, permute(meanDenResp,[3,1,2])});
cl = squeeze(ci(1,:,:,:));
cu = squeeze(ci(2,:,:,:));

MakeFigure;
PlotAsymmetricErrorPatch(updateRate, mu, cl,cu, corder);
legend(legendStr);
xlabel('update rate (Hz)');
ylabel('denominator LNLN response (arb. units)');
ConfAxis(16);
title('Denominator LNLN');
axis('square');

%% Plot the numerator and denominator LNLN responses relative to uncorrelated noise

meanUn = mean(meanNumResp(:,3,:),3);
relResp = (meanNumResp - meanUn) ./ meanUn;
mu = nanmean(relResp,3);
ci = bootci(params.nboot, {@nanmean, permute(relResp,[3,1,2])});
cl = squeeze(ci(1,:,:,:));
cu = squeeze(ci(2,:,:,:));

MakeFigure;
PlotAsymmetricErrorPatch(updateRate, mu, cl,cu, corder);
legend(legendStr);
xlabel('update rate (Hz)');
ylabel('numerator LNLN response relative to uncorrelated');
a = ceil(max(abs(relResp(:))));
ylim([-a,a]);
ConfAxis(16);
title('Numerator LNLN');
axis('square');

meanUn = mean(meanDenResp(:,3,:),3);
relResp = (meanDenResp - meanUn) ./ meanUn;
mu = nanmean(relResp,3);
ci = bootci(params.nboot, {@nanmean, permute(relResp,[3,1,2])});
cl = squeeze(ci(1,:,:,:));
cu = squeeze(ci(2,:,:,:));

MakeFigure;
PlotAsymmetricErrorPatch(updateRate, mu, cl,cu, corder);
legend(legendStr);
xlabel('update rate (Hz)');
ylabel('denominator LNLN response relative to uncorrelated');
a = ceil(max(abs(relResp(:))));
ylim([-a,a]);
ConfAxis(16);
title('Denominator LNLN');
axis('square');

%% Plot the numerator and denominator LNLN DSIs

dsiNum = (meanNumResp(:,[1,3],:)-meanNumResp(:,[2,4],:)) ./ (meanNumResp(:,[1,3],:)+meanNumResp(:,[2,4],:));
mu = nanmean(dsiNum,3);
ci = bootci(params.nboot, {@nanmean, permute(dsiNum,[3,1,2])});
cl = squeeze(ci(1,:,:,:));
cu = squeeze(ci(2,:,:,:));

MakeFigure;
PlotAsymmetricErrorPatch(updateRate, mu, cl,cu, corder);
legend({'+','-'});
xlabel('update rate (Hz)');
ylabel('numerator LNLN DSI');
ConfAxis(16);
axis('square');
title('Numerator LNLN');
ylim([-1,1]);

dsiDen = (meanDenResp(:,[1,3],:)-meanDenResp(:,[2,4],:)) ./ (meanDenResp(:,[1,3],:)+meanDenResp(:,[2,4],:));
mu = nanmean(dsiDen,3);
ci = bootci(params.nboot, {@nanmean, permute(dsiDen,[3,1,2])});
cl = squeeze(ci(1,:,:,:));
cu = squeeze(ci(2,:,:,:));

MakeFigure;
PlotAsymmetricErrorPatch(updateRate, mu, cl,cu, corder);
legend({'+','-'});
xlabel('update rate (Hz)');
ylabel('denominator LNLN DSI');
ConfAxis(16);
axis('square');
title('Denominator LNLN');
ylim([-1,1]);

%% Plot all DSIs together

dsiAll = cat(2, dsi(:,1:2,:), dsiNum(:,1:2,:), dsiDen(:,1:2,:));
mu = nanmean(dsiAll,3);
ci = bootci(params.nboot, {@nanmean, permute(dsiAll,[3,1,2])});
cl = squeeze(ci(1,:,:,:));
cu = squeeze(ci(2,:,:,:));

MakeFigure;
PlotAsymmetricErrorPatch(updateRate, mu, cl,cu, corder2(1:6,:));
legend({'full model +','full model -','numerator +','numerator -','denominator +','denominator -'});
xlabel('update rate (Hz)');
ylabel('DSI');
ylim([-1,1]);
ConfAxis(16);
axis('square');

end