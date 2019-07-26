function ThreeInputModelLinearReceptiveField(config, varargin)

%% Either load data from cache or run simulation

% Search for existing cache files
dataPath = dir(fullfile(config.rootDataPath, 'linearReceptiveFieldData_*.mat'));

% Set list of variable names to save or load
saveVarList = {'params','filters','noiseParam','kernel','tau','rho'};

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
    
    % Kernel parameters
    extentForwards = 60;
    extentBackwards = 0;
    shiftX = 50;
    shiftDx = params.dx;
    
    % Define sampling information
    idx = params.t>params.tAv;
    isSampled = rem((1:nnz(idx))',round(1/(params.dt * noiseParam.updateRate))) == 0;
    sampleIdxs = (1:nnz(isSampled))'-1;
    validSampleIdxs = sampleIdxs( (sampleIdxs > extentForwards) & (sampleIdxs < length(sampleIdxs)-extentBackwards+1));
    numX = length(params.x);
    shiftSpacing = floor(shiftDx/params.dx);
    numShift = floor(shiftX/shiftDx);
    
    
    % Define plotting vectors
    tau = (-extentBackwards:extentForwards)' / noiseParam.updateRate;
    rho = (-numShift:shiftSpacing:numShift)' * params.dx;
    
    %% Run the simulation
    
    % Allocate a container
    kernel = nan(extentBackwards + extentForwards + 1, 2*numShift+1,params.numRep);
    
    % Iterate over realizations
    parfor indR = 1:params.numRep
        tic;
        
        % Make the stimulus
        [ stimArray ] = UncorrelatedBinaryNoise(params, noiseParam);
        
        % Compute the response
        [ meanResp, voltageResp, calciumResp ] = ComputeThreeInputModelResponse(stimArray, params, filters);
        
        % Sample the stimulus and response
        stimAll = stimArray(isSampled, :);
        respAll = calciumResp(isSampled,:);
        
        % Allocate a container to store the local kernel estimate
        locKernel = nan(extentBackwards + extentForwards + 1, 2*numShift+1, numX);
        
        % Iterate over spatial shifts
        for indS = -numShift:numShift
            
            % Shift the stimulus
            stimShift = circshift(stimAll, -indS*shiftSpacing, 2);
            
            % Iterate over positions in space
            for indX = 1:numX
                
                % Grab the stimulus and response
                stim = stimShift(:,indX);
                validSampledResp = respAll(validSampleIdxs,indX);
                validSampledResp = validSampledResp - mean(validSampledResp);
                stimMatrix = calcSparseToeplitz(stim,sampleIdxs,extentForwards,extentBackwards);
                stimMatrix = stimMatrix-mean(stimMatrix);
                
                % Compute the kernel
                if params.olsKernel
                    % If desired, extract using OLS
                    locKernel(:,indS+numShift+1,indX) = pinv(stimMatrix) * validSampledResp;
                else
                    % Otherwise, use reverse correlation
                    locKernel(:,indS+numShift+1,indX) = noiseParam.updateRate .* mean(stimMatrix .* validSampledResp,1) ./ var(stim,0,1);
                end
            end
        end
        
        % Average over space
        kernel(:,:,indR) = mean(locKernel,3);
        
        % Print a status update
        fprintf('Realization %d of %d: %f seconds.\n', indR, params.numRep, toc);
    end
    
    % Save the data
    tic;
    savePath = fullfile(config.rootDataPath, sprintf('linearReceptiveFieldData_%s.mat', datestr(datetime('now'), 'yyyymmdd_HHMMSS')));
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

load('utils/blueRedColorMap.mat', 'cmpRed','cmpBlueRed');

%% Compute kernel statistics

% Average the kernel over realizations
meanKernel = mean(kernel,3);

% % Compute CI for the kernel by bootstrapping over realizations of the noise
% ciKernel = bootci(params.nboot, {@mean, permute(kernel, [3,1,2])});
% ciLowerKernel = squeeze(ciKernel(1,:,:,:));
% ciUpperKernel = squeeze(ciKernel(2,:,:,:));

%% Show the linear RF

% Plot full linear RF
MakeFigure;
imagesc(rho, tau, meanKernel);
colormap(cmpBlueRed);
kMax = max(abs(meanKernel(:)));
caxis([-kMax,kMax]);
xlabel('relative azimuthal location');
ylabel('time in past (s)');
yticks(0:0.25:1);
axis('square');
cbar = colorbar;
ylabel(cbar, 'filter strength (arb. units/s/c)');
ConfAxis(16);

%% Plot an example timeseries

% Make the stimulus
[ stimArray ] = UncorrelatedBinaryNoise(params, noiseParam);

% Compute the response
[ meanResp, voltageResp, calciumResp ] = ComputeThreeInputModelResponse(stimArray, params, filters);

% Find spatial central point
[~,spatialInd] = min(abs(params.x - 180));

MakeFigure;
subplot(3,1,1);
plot(params.t, stimArray(:,spatialInd), '-','linewidth', 2);
xlabel('time (s)');
ylabel('contrast');
ylim([-1.1,1.1]);
ConfAxis(16);

subplot(3,1,2);
plot(params.t, voltageResp(:,spatialInd), '-','linewidth', 2);
xlabel('time (s)');
ylabel('V_{mem} (mV)');
ylim([-15 15]);
yticks([-15 0 15]);
ConfAxis(16);

subplot(3,1,3);
plot(params.t, calciumResp(:,spatialInd), '-','linewidth', 2);
xlabel('time (s)');
ylabel('response (arb. units)');
ConfAxis(16);

end
