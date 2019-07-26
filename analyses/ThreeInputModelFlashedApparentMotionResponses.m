function ThreeInputModelFlashedApparentMotionResponses()

% Set overall parameters
[ params ] = SetModelParameters('tInt', 1/2, 'tAv', 0, 'tOn', 2, 'dx', 0.1);
[ filters ] = MakeModelFilters(params);

% Set stimulus parameters
barParam.mlum = 0;
barParam.c = 1;
barParam.duration = 0.4;
barParam.barWidth = 4.5;

%% Run the simulation
tic;

% Set bar center location
barCenter = 180;

% Present all components and two-bar flashed apparent motion pairs
[ stimArray ] = TwoBarFlashedApparentMotion(params, barParam, barCenter);

% Compute responses
[ meanResp, voltageResp, calciumResp ] = ComputeThreeInputModelResponse(stimArray, params, filters);

fprintf('Completed simulation in %f seconds.\n', toc);

%% Set options for plotting

barPairLegendStr = {'++PD','++ND','--PD','--ND','+-PD','+-ND','-+PD','-+ND'};
% corder = lines(4);
corder = lines(1);

%% Plot kymographs of stimuli

MakeFigure;
for ind = 1:8
    subplot(2,4,ind);
    imagesc(params.t, params.x-barCenter, stimArray(:,:,ind)');
    colormap([0,0,0;0.5,0.5,0.5;1,1,1]);
    caxis([-1,1]);
    xlim([-0.25, 1]);
    ylim([-20 20]);
    % title(barPairLegendStr{ind});
    axis('xy','square');
    xlabel('time (s)');
    ylabel('relative azimuthal location (\circ)');
end

MakeFigure;
for ind = 1:8
    subplot(2,4,ind);
    imagesc(params.t, params.x-barCenter, stimArray(:,:,8+ind)');
    colormap([0,0,0;0.5,0.5,0.5;1,1,1]);
    caxis([-1,1]);
    xlim([-0.25, 1]);
    ylim([-20 20]);
    title(barPairLegendStr{ind});
    axis('xy','square');
    xlabel('time (s)');
    ylabel('relative azimuthal location (\circ)');
end

%% Plot the results

% Find spatial central point
[~,spatialInd] = min(abs(params.x - 180));

% Extract the desired response
calciumRespSel = squeeze(calciumResp(:,spatialInd,:));

linearPrediction = nan(length(params.t),8);
linearPrediction(:,1) = sum(calciumRespSel(:,[1,2]),2);
linearPrediction(:,2) = sum(calciumRespSel(:,[3,4]),2);
linearPrediction(:,3) = sum(calciumRespSel(:,[5,6]),2);
linearPrediction(:,4) = sum(calciumRespSel(:,[7,8]),2);
linearPrediction(:,5) = sum(calciumRespSel(:,[1,6]),2);
linearPrediction(:,6) = sum(calciumRespSel(:,[3,8]),2);
linearPrediction(:,7) = sum(calciumRespSel(:,[2,5]),2);
linearPrediction(:,8) = sum(calciumRespSel(:,[4,7]),2);

MakeFigure;
plot(params.t, calciumRespSel(:,[9:10])-linearPrediction(:,1:2), 'linewidth', 2);
legend(barPairLegendStr(1:2));
xlabel('time (s)');
ylabel('nonlinear response component (arb. units)');
title(sprintf('%0.2f\\circ bars, %0.2f seconds',barParam.barWidth,barParam.duration));
ConfAxis(16);

MakeFigure;
hold on;
bar(1:8, mean(calciumRespSel(:,9:16)), 'facecolor','flat','cdata', corder, 'EdgeColor','none');
bar(1:8, mean(linearPrediction), 'facecolor','none','LineWidth',2);
xlim([0 9]);
xticks(1:8);
xticklabels(barPairLegendStr);
ylabel('response (arb. units)');
title(sprintf('%0.2f\\circ bars, %0.2f seconds',barParam.barWidth,barParam.duration));
ConfAxis(16);

MakeFigure;
hold on;
bar(1:8, mean(calciumRespSel(:,9:16)) ./ mean(linearPrediction), 'facecolor','flat','cdata', corder, 'EdgeColor','none');
plot([0 9], [1,1], '--k','linewidth',2);
xlim([0 9]);
xticks(1:8);
xticklabels(barPairLegendStr);
ylabel('response relative to linear prediction');
title(sprintf('%0.2f\\circ bars, %0.2f seconds',barParam.barWidth,barParam.duration));
ConfAxis(16);

end
