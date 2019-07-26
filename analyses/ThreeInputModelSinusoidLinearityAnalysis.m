function ThreeInputModelSinusoidLinearityAnalysis(config, varargin)

% Set overall parameters
%     [ params ] = SetModelParameters('tInt', 0, varargin{:});
[ params ] = SetModelParameters('tInt', 0, 'tOn', 3, 'dx', 1);

% Make the filters
[ filters ] = MakeModelFilters(params);

% Contrast
% c = 1/2;
c = 1;

% Temporal frequency
tf = 1;

% Spatial wavelength
lambda = 25;
% lambda = 45;

% Set line colors
% corder = [0.9153    0.2816    0.2878;0.3467    0.5360    0.6907];
corder = lines(2);

% Find spatial central point
[~,spatialInd] = min(abs(params.x - 180));

%% Run the simulation
tic;

% Make the stimuli
[ stimArrayM ] = MonocomponentSinusoidalGratings(params, c, tf, 1/lambda);
[ stimArrayC ] = CompositeSinusoidalGratingsForLinearityAnalysis(params, c, tf, 1/lambda);

% Compute the responses
[ meanRespM, voltageRespM, calciumRespM ] = ComputeThreeInputModelResponse(stimArrayM, params, filters);
[ meanRespC, voltageRespC, calciumRespC ] = ComputeThreeInputModelResponse(stimArrayC, params, filters);

fprintf('Completed simulation in %f seconds.\n', toc);

%% Compute the linear predictions

t = params.t(params.averagingMask);
respReal = squeeze(voltageRespM(params.averagingMask,spatialInd,:));
respPred = squeeze(sum(voltageRespC(params.averagingMask,spatialInd,:,:),4))/4;

%% Plot the results

minMaxV = [-40,40];

% Plot the PD predictions
MakeFigure;
plot(t, respReal(:,1), '-', 'linewidth', 2, 'color', corder(1,:));
hold on;
plot(t, respPred(:,1), '-', 'linewidth', 2, 'color', [0.5,0.5,0.5]);
axis('off');
plot([0, 0, .5], [10, 0, 0], '-k', 'linewidth', 2);
h = text(-0.05, 5, '10 mV', 'horiz','right','vert','middle', 'fontsize', 16);
set(h, 'rotation', 90);
text(0.15, -0.05, '500 ms', 'horiz','center','vert','top', 'fontsize', 16);
ylim(minMaxV);
legend({'Model PD response','Linear prediction of PD response'},'location','northwest');
legend boxoff;
ConfAxis(16);
title(sprintf('%0.2f Hz, %0.2f\\circ, contrast %0.2f', tf, lambda, c));

% Plot the ND predictions
MakeFigure;
plot(t, respReal(:,2), '-', 'linewidth', 2, 'color', corder(2,:));
hold on;
plot(t, respPred(:,2), '-', 'linewidth', 2, 'color', [0.5,0.5,0.5]);
axis('off');
plot([0, 0, .5], [10, 0, 0], '-k', 'linewidth', 2);
h = text(-0.05, 5, '10 mV', 'horiz','right','vert','middle', 'fontsize', 16);
set(h, 'rotation', 90);
text(0.15, -0.05, '500 ms', 'horiz','center','vert','top', 'fontsize', 16);
ylim(minMaxV);
legend({'Model ND response','Linear prediction of ND response'},'location','northwest');
legend boxoff;
ConfAxis(16);
title(sprintf('%0.2f Hz, %0.2f\\circ, contrast %0.2f', tf, lambda, c));

%% Scatter predicted and actual responses against one another

r2 = 1-sum((respPred-respReal).^2,1)./sum((respReal-mean(respReal,1)).^2,1);

MakeFigure;
scatter(respPred(:,1), respReal(:,1), 20, corder(1,:), 'filled');
hold on;
scatter(respPred(:,2), respReal(:,2), 20, corder(2,:), 'filled');
plot(minMaxV, minMaxV, '--k','linewidth',2);
axis('equal');
xlim(minMaxV);
ylim(minMaxV);
xlabel('linear prediction (mV)');
ylabel('model response (mV)');
legend({sprintf('PD (R^2 = %0.2f)', r2(1)),sprintf('ND (R^2 = %0.2f)', r2(2))}, 'location','northwest');
ConfAxis(16);
title(sprintf('%0.2f Hz, %0.2f\\circ, contrast %0.2f', tf, lambda, c));

end