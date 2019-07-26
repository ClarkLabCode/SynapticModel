function ThreeInputModelMovingEdgeResponses(config, varargin)

%% Set parameters

% [ params ] = SetModelParameters('tInt', 1, 'tAv', 0, varargin{:});
[ params ] = SetModelParameters('tInt', 2, 'tAv', 0, 'tOn', 12);
[ filters ] = MakeModelFilters(params);

%%
% Velocities
v = (30)';

% Stimulus names
legendStr = {'PD ON','PD OFF','ND ON','ND OFF'};

% corder = [1,0,0; 0,0,1; 1,0.64,0.5; 0,1,0];
corder = lines(4);

% barParam.barWidth = 5;
% barParam.barPeriod = 45;
barParam.mlum = 0;
barParam.c = 1;

%% Run the simulation
tic;

[ stimArray ] = MovingEdges(params, barParam, v);
[ meanResp, voltageResp, calciumResp ] = ComputeThreeInputModelResponse(stimArray, params, filters);

fprintf('Completed simulation in %f seconds.\n', toc);

%% Plot voltage response timeseries

% Find spatial central point
[~,spatialInd] = min(abs(params.x - 180));

sp = 5;
MakeFigure;
hold on;
for ind = 1:4
    plot(params.t(params.mask), squeeze(voltageResp(params.mask,spatialInd,ind))-(ind-1)*sp, 'linewidth', 2, 'color', corder(ind,:));
end
xlim([0 10]);
axis('off');
plot([0, 0, .5], [10, 0, 0], '-k', 'linewidth', 2);
h = text(-0.2, 8, '10 mV', 'horiz','right','vert','middle', 'fontsize', 16);
set(h, 'rotation', 90);
text(0.15, 0, '2 s', 'horiz','center','vert','top', 'fontsize', 16);
legend(legendStr);
legend('boxoff');
ConfAxis(16);

%% Plot calcium response timeseries

sp = 50;
MakeFigure;
hold on;
for ind = 1:4
    plot(params.t, squeeze(calciumResp(:,spatialInd,ind))-(ind-1)*sp, 'linewidth', 2, 'color', corder(ind,:));
end
xlim([0 10]);
axis('off');
plot([0, 0, 2], [150, 50, 50], '-k', 'linewidth', 2);
h = text(-0.2, 140, '100 arb. units', 'horiz','right','vert','middle', 'fontsize', 16);
set(h, 'rotation', 90);
text(0.15, 50, '2 s', 'horiz','center','vert','top', 'fontsize', 16);
legend(legendStr);
legend('boxoff');
ConfAxis(16);

end