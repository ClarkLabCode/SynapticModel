function ThreeInputModelLocalFlashedBarResponse()

%% Set parameters

% [ params ] = SetModelParameters('tInt', 0, varargin{:});
[ params ] = SetModelParameters('tInt', 1/2, 'dt', 1e-3, 'dx', 0.1);
[ filters ] = MakeModelFilters(params);

barParam.barWidth = 2;
barParam.mlum = 0;
barParam.c = 1;

% barDuration = [0.02; 0.04; 0.08; 0.160];
barDuration = 0.160;

load('utils/blueRedColorMap.mat', 'cmpRed','cmpBlueRed');

%% Compute the response

tic;

[ stimArray ] = LocalFlashedBar(params, barParam, barDuration);

[ meanResp, voltageResp, calciumResp ] = ComputeThreeInputModelResponse(stimArray, params, filters);

fprintf('Completed simulation in %f seconds.\n', toc);

%% Plot the results

[~,centerInd] = min(abs(params.x-mean(params.x)));

xExtent = floor(15 / params.dx);
x = (-xExtent:xExtent)' * params.dx;

% Select desired spatial domain
voltageRespSel = fliplr(voltageResp(:,(-xExtent:xExtent)+centerInd,:));
stimArraySel = stimArray(:, (-xExtent:xExtent)+centerInd,:);

MakeFigure;
imagesc(x, params.t, stimArraySel(:,:,1));
ylim([-0.1,1])
cbar = colorbar;
caxis([-1 1])
colormap([0,0,0; 0.5,0.5,0.5; 1,1,1]);
axis('square');
xlabel('relative azimuthal location(\circ)');
ylabel('time (s)');
ylabel(cbar, 'contrast');
ConfAxis(16);
title(sprintf('%d ms, %d{\\circ} flashed white bar on gray background', barDuration*1000, barParam.barWidth));


MakeFigure;
imagesc(x, params.t, voltageRespSel(:,:,1));
hold on;
contour(x, params.t, voltageRespSel(:,:,1), -15:5:15, 'EdgeColor','k','LineWidth',2);
ylim([-0.1,1]);
plot(x, zeros(length(x),1), '--k','linewidth',2);
plot(x, barDuration*ones(length(x),1), '--k','linewidth',2);
cbar = colorbar;
caxis([-15 15])
colormap(cmpBlueRed)
axis('square');
xlabel('relative azimuthal location(\circ)');
ylabel('time (s)');
ylabel(cbar, 'membrane voltage (mV)');
ConfAxis(16);
title(sprintf('%d ms, %d{\\circ} flashed white bar on gray background', barDuration*1000, barParam.barWidth));

MakeFigure;
imagesc(x, params.t, stimArraySel(:,:,2));
ylim([-0.1,1])
cbar = colorbar;
caxis([-1 1])
colormap([0,0,0; 0.5,0.5,0.5; 1,1,1]);
axis('square');
xlabel('relative azimuthal location(\circ)');
ylabel('time (s)');
ylabel(cbar, 'contrast');
ConfAxis(16);
title(sprintf('%d ms, %d{\\circ} flashed black bar on gray background', barDuration*1000, barParam.barWidth));

MakeFigure;
imagesc(x, params.t, voltageRespSel(:,:,2));
hold on;
contour(x, params.t, voltageRespSel(:,:,2), -15:5:15, 'EdgeColor','k','LineWidth',2);
ylim([-0.1,1]);
plot(x, zeros(length(x),1), '--k','linewidth',2);
plot(x, barDuration*ones(length(x),1), '--k','linewidth',2);
cbar = colorbar;
caxis([-15 15])
colormap(cmpBlueRed)
axis('square');
xlabel('relative azimuthal location(\circ)');
ylabel('time (s)');
ylabel(cbar, 'membrane voltage (mV)');
ConfAxis(16);
title(sprintf('%d ms, %d{\\circ} flashed black bar on gray background', barDuration*1000, barParam.barWidth));

end
