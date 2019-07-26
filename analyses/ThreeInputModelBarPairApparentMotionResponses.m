function ThreeInputModelBarPairApparentMotionResponses(config, varargin)

%% Either load data from cache or run simulation

% Search for existing cache files
dataPath = dir(fullfile(config.rootDataPath, 'barPairApparentMotionData_*.mat'));

% Set list of variable names to save or load
saveVarList = {'params','filters','barParam','singleBarStimArray','barPairStimArray','singleBarVoltageResp','singleBarCalciumResp','barPairVoltageResp','barPairCalciumResp'};

% Check if data regeneration flag is set or if no cache exists
if config.regenerateData || isempty(dataPath)
    fprintf('Simulation run flag set or no cache found. Running simulation...\n');
    tic;

    % Set overall parameters
    [ params ] = SetModelParameters('tInt', 1, 'tAv', 0, 'tOn', 1, 'dx', 0.1, varargin{:});
    [ filters ] = MakeModelFilters(params);

    % Set stimulus parameters
    barParam.mlum = 0;
    barParam.c = 1;
    barParam.delay = 0.15;
    barParam.barWidth = 5;
    barParam.barOffset = 5;
    barParam.barPeriod = 45;

    % Make stimuli

    % Full leading bar
    [ singleBarStimArray ] = SingleBars(params, barParam);

    % Bar pair
    [ barPairStimArray ] = BarPairs(params, barParam);

    % Compute single-bar responses
    [ singleBarMeanResp, singleBarVoltageResp, singleBarCalciumResp ] = ComputeThreeInputModelResponse(singleBarStimArray, params, filters);

    % Compute bar pair responses
    [ barPairMeanResp, barPairVoltageResp, barPairCalciumResp ] = ComputeThreeInputModelResponse(barPairStimArray, params, filters);

    % Print a status update
    fprintf('Completed simulation in %f seconds\n', toc);

    % Save the data
    tic;
    savePath = fullfile(config.rootDataPath, sprintf('barPairApparentMotionData_%s.mat', datestr(datetime('now'), 'yyyymmdd_HHMMSS')));
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

singleBarLegendStr = {'+','-'};
barPairLegendStr = {'++PD','++ND','--PD','--ND','+-PD','+-ND','-+PD','-+ND'};

load('utils/blueRedColorMap.mat');

%% Average responses over time

% Single bar average
singleBarMeanRespOverTime = squeeze(nanmean(singleBarCalciumResp(params.averagingMask,:,:),1));

% Make a temporal indexing vector
timeMask = params.averagingMask & params.t > barParam.delay;

%% Align bar pair responses

detectorLoc = params.x;

% Align to peak single bar response
x = max(singleBarMeanRespOverTime, [], 2);
[~,barIdxSingleBar] = findpeaks(x, 'minPeakHeight', max(x)/2);

% Average desired responses over space
barPairSingleBarAlignedMaskResp = squeeze(nanmean(barPairCalciumResp(:,barIdxSingleBar,:),2));

% Compute the average responses
barPairSingleBarAlignedResp = squeeze(nanmean(barPairSingleBarAlignedMaskResp(timeMask,:),1));

%% Show kymograph of single-bar stimuli

% Plot at full resolution
for ind = 1:size(singleBarStimArray,3)
    MakeFigure;
    imagesc(params.t, params.x, singleBarStimArray(:,:,ind)');
    axis('xy','square','tight');
    xlabel('time (s)');
    ylabel('spatial position (\circ)');
    ylim([0 params.xExtent]);
    cbar = colorbar;
    ylabel(cbar, 'input contrast');
    cbar.Ticks = [-1,0,1];
    colormap([0,0,0;1/2,1/2,1/2;1,1,1]);
    % colormap(cmpBlueRed);
    ConfAxis(16);
    caxis([-barParam.c,barParam.c]+barParam.mlum);
    title(singleBarLegendStr{ind});
    ylim([0 2*barParam.barPeriod]);
end

%% Show kymograph of single-bar voltage responses

for ind = 1:size(singleBarVoltageResp,3)
    MakeFigure;
    imagesc(params.t, params.x, singleBarVoltageResp(:,:,ind)');
    axis('xy','square','tight');
    xlabel('time (s)');
    ylabel('spatial position (\circ)');
    cbar = colorbar;
    ylabel(cbar, 'membrane voltage (mV)');
    caxis([-20 20])
    colormap(cmpBlueRed)
    ConfAxis(16);
    hold on;
    plot(zeros(length(params.x),1), params.x, '--k', 'linewidth',2);
    plot(params.tOn*ones(length(params.x),1), params.x, '--k', 'linewidth',2);
    title(singleBarLegendStr{ind});
    ylim([0 2*barParam.barPeriod]);
end

%% Show kymograph of single-bar calcium responses

for ind = 1:size(singleBarVoltageResp,3)
    MakeFigure;
    imagesc(params.t, params.x, singleBarCalciumResp(:,:,ind)');
    axis('xy','square','tight');
    xlabel('time (s)');
    ylabel('spatial position (\circ)');
    cbar = colorbar;
    ylabel(cbar, 'response (arb. units)');
    colormap(cmpRed)
    ConfAxis(16);
    hold on;
    plot(zeros(length(params.x),1), params.x, '--k', 'linewidth',2);
    plot(params.tOn*ones(length(params.x),1), params.x, '--k', 'linewidth',2);
    title(singleBarLegendStr{ind});
    ylim([0 2*barParam.barPeriod]);
end

%% Plot single bar responses averaged over time

MakeFigure;
plot(detectorLoc, singleBarMeanRespOverTime, 'linewidth', 2);
legend(singleBarLegendStr);
axis('square');
xlabel('spatial position (\circ)');
ylabel('response (arb. units)');
title('single bars');
ConfAxis(16);
% ylim([0 1.2]);

%% Plot time-averaged bar pair responses over space

% barPaiMeanRespOverTime = squeeze(nanmean(barPairCalciumResp(timeMask,:,:),1));

% barPairMeanRespOverTime = squeeze(nanmean(barPairCalciumResp(timeMask,:,:) / normFactor,1));
barPairMeanRespOverTime = squeeze(nanmean(barPairCalciumResp(timeMask,:,:),1));
legendStr = {'phi', 'reverse-phi'};
for ind = 1:4:size(barPairMeanRespOverTime,2)
    MakeFigure;
    hold on;
    plot(detectorLoc, barPairMeanRespOverTime(:, ind:ind+3), 'linewidth', 2);
    legend(barPairLegendStr(ind:ind+3));
    title(legendStr{floor(ind/4)+1});
    axis('square');
    xlabel('spatial location (\circ)');
    ylabel('response (arb. units)');
    ConfAxis(16);
    %     ylim([0 1.2]);
end

%% Plot timeseries of bar pair responses aligned to peak single bar

legendStr = {'phi','reverse-phi'};
for ind = 1:4:size(barPairSingleBarAlignedMaskResp,2)
    MakeFigure;
    hold on;
    plot(params.t, barPairSingleBarAlignedMaskResp(:, ind:ind+3), 'linewidth', 2);
    plot([0,0], [0,200], '--k', 'linewidth', 2);
    plot([barParam.delay,barParam.delay], [0,200], '--k', 'linewidth', 2);
    plot([params.tOn,params.tOn], [0,200], '--k', 'linewidth', 2);
    legend(barPairLegendStr(ind:ind+3));
    title(strcat(legendStr{floor(ind/4)+1}, ', aligned to peak single bar response'));
    axis('square');
    xlabel('time (s)');
    ylabel('response (arb. units)');
    ConfAxis(16);
    %     ylim([0 1.2]);
end


%% Make a bar graph of bar pair responses

MakeFigure;
bar(barPairSingleBarAlignedResp / max(barPairSingleBarAlignedResp([1,3])));
% bar(barPairResp);
xticks(1:8);
xlim([0 9]);
xticklabels(barPairLegendStr);
ylabel('response (arb. units)');
xlabel('bar pairing');
ylim([0 1.2]);
axis('square');
ConfAxis(16);

end
