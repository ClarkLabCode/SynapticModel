function ConvertNaturalScenesToContrast(config, params)
%ConvertNaturalScenesToContrast: Converts natural scene database to
%contrast and saves results to a data file.

%% Load natural scenes from file

tic;
files = dir(fullfile(config.sceneSourcePath, 'imageData','*.mat'));
scenes = [];
for ind = 1:size(files,1)
    fpath = fullfile(files(ind).folder, files(ind).name);
    data = load(fpath);
    scenes = cat(3, scenes, data.projection);
    fprintf('\tLoaded file: %s\n', fpath);
end
fprintf('Loaded natural scene data in %f seconds.\n', toc);

%% Make blurring filters

% Photoreceptor kernel sd
filtStd = params.fwhmBlur/(2*sqrt(2*log(2)));

% Mean estimation kernel sd
filtStdContrast = params.fwhmAverage/(2*sqrt(2*log(2)));

%% Convert scenes to contrast
tic;

% Get the spatial resolution
xRes = 360/size(scenes,2);

% Make x and y vectors
x = (0:xRes:360-xRes)';
y = (0:xRes:xRes*size(scenes,1)-xRes)';

% Define photoreceptor blurring kernels
xFilt = ifftshift(normpdf(x,180,filtStd));
yFilt = ifftshift(normpdf(y,(y(end)+xRes)/2,filtStd));
xyFiltMat = yFilt*xFilt';
xyFiltTensor = repmat(xyFiltMat,[1 1 size(scenes,3)]);

% Define contrast filter kernels
xFiltContrast = ifftshift(normpdf(x,180,filtStdContrast));
yFiltContrast = ifftshift(normpdf(y,(y(end)+xRes)/2,filtStdContrast));
xyFiltMatContrast = yFiltContrast*xFiltContrast';
xyFiltTensorContrast = repmat(xyFiltMatContrast,[1 1 size(scenes,3)]);

% Compute the Fourier transforms of the filters and scenes
fftScenes = fft2(scenes);
fftFilt = fft2(xyFiltTensor);
fftFiltContrast = fft2(xyFiltTensorContrast);

% Filter the scenes
filteredScenes = real(ifft2(fftScenes.*fftFilt));
filteredScenesContrast = real(ifft2(fftScenes.*fftFiltContrast));

% Compute the contrast
contrastScenes = (filteredScenes-filteredScenesContrast)./filteredScenesContrast;

% Print a status update
fprintf('Converted natural scenes to contrast in %f seconds.\n', toc);

%% Save natural scene data to file
tic;

% Form the file name
savePath = sprintf('naturalSceneData_%s.mat', datestr(datetime('now'), 'yyyymmdd_HHMMSS'));

% Save the contrast scene data
save(fullfile(config.rootDataPath, savePath), 'contrastScenes', '-v7.3');

% Print a status update
fprintf('Saved contrast scenes in %f seconds.\n', toc);

end

