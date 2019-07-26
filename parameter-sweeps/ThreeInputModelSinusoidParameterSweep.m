function ThreeInputModelSinusoidParameterSweep(gExc, gInh, varargin)



%% Set the grids

if nargin < 1 || isempty(gExc)
    gExc = (0.05:0.05:1)';
end
if nargin < 2 || isempty(gInh)
    gInh = (0.05:0.05:1)';
end

%% Set overall parameters

[ paramsDefault ] = SetModelParameters(varargin{:});
[ filters ] = MakeModelFilters(paramsDefault);

%% Set stimulus parameters

% Contrast
c = 1/2;

% Spatial frequency
sf = 1/45;

% Temporal frequency
tf = 1;

legendStr = {'PD','ND','PD+ND','PD+OD'};

%% Set options for plotting

% Set colormap
load('utils/blueRedColorMap.mat', 'cmpRed','cmpBlueRed');

% Set options for overlaid contours
contourOpts = {-1:0.1:1,'EdgeColor', 'k', 'linewidth', 2};

%% Make stimuli

% Make the monocomponent sinusoid stimulus array
[ monoStimArray ] = MonocomponentSinusoidalGratings(paramsDefault, c, tf, sf);

% Make the composite sinusoid stimulus array
[ compStimArray ] = CompositeSinusoidalGratings(paramsDefault, c, tf, sf);

%% Run the parameter sweep

nExc = length(gExc);
nInh = length(gInh);

monoMeanResp = nan(nInh, nExc, 2);
compMeanResp = nan(nInh, nExc, 2);

for indE = 1:nExc
    
    parfor indI = 1:nInh
        
        tic;
        
        % Grab the parameter values
        params = paramsDefault;
        params.g1 = gInh(indI);
        params.g2 = gExc(indE);
        params.g3 = gInh(indI);
        
        % Compute monocomponent sinusoid responses
        [ monoMeanResp(indI, indE, :) ] = ComputeThreeInputModelResponse(monoStimArray, params, filters);
        
        % Compute composite sinusoid responses
        [ compMeanResp(indI, indE, :) ] = ComputeThreeInputModelResponse(compStimArray, params, filters);
        
        fprintf('gExc %d of %d, gInh %d of %d: %f s\n', indE, nExc, indI, nInh, toc);
    end
end

% Combine everything together
meanResp = cat(3, monoMeanResp, compMeanResp);

%% Plot raw responses

for ind = 1:size(meanResp,3)
    MakeFigure;
    
    resp = meanResp(:,:,ind) / max(meanResp(:));
    
    imagesc(gExc, gInh, resp);
    hold on;
    contour(gExc, gInh, resp, contourOpts{:});
    axis('xy','square','tight');
    caxis([0,1]);
    cbar = colorbar;
    
    ylabel(cbar, 'response (arb. units)');
    colormap(cmpRed);
    title(legendStr{ind});
    xlabel('g_{exc} / g_{leak}');
    ylabel('g_{inh} / g_{leak}');
    ConfAxis(16);
end

%% Plot response indices relative to PD

% DSI
dsi = (meanResp(:,:,1) - meanResp(:,:,2)) ./ (meanResp(:,:,1)+meanResp(:,:,2));
MakeFigure;
imagesc(gExc, gInh, dsi);
hold on;
contour(gExc, gInh, dsi, contourOpts{:});
axis('xy','square','tight');
caxis([-1,1]);
cbar = colorbar;
ylabel(cbar, 'DSI');
colormap(cmpBlueRed);
title('DSI');
xlabel('g_{exc} / g_{leak}');
ylabel('g_{inh} / g_{leak}');
ConfAxis(16);

% Opponency index
dsi = (meanResp(:,:,3) - meanResp(:,:,1)) ./ (meanResp(:,:,3)+meanResp(:,:,1));
MakeFigure;
imagesc(gExc, gInh, dsi);
hold on;
contour(gExc, gInh, dsi, contourOpts{:});
axis('xy','square','tight');
caxis([-1,1]);
cbar = colorbar;
ylabel(cbar, 'opponency index');
colormap(cmpBlueRed);
title('opponency index');
xlabel('g_{exc} / g_{leak}');
ylabel('g_{inh} / g_{leak}');
ConfAxis(16);

% OD suppression index
dsi = (meanResp(:,:,4) - meanResp(:,:,1)) ./ (meanResp(:,:,4)+meanResp(:,:,1));
MakeFigure;
imagesc(gExc, gInh, dsi);
hold on;
contour(gExc, gInh, dsi, contourOpts{:});
axis('xy','square','tight');
caxis([-1,1]);
cbar = colorbar;
ylabel(cbar, 'OD suppression index');
colormap(cmpBlueRed);
title('OD suppression index');
xlabel('g_{exc} / g_{leak}');
ylabel('g_{inh} / g_{leak}');
ConfAxis(16);

end
