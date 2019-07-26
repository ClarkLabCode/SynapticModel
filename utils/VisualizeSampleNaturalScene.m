
%%

files = dir(fullfile(config.sceneSourcePath, 'imageData','*.mat'));
scenes = [];
ind = randi([1, 421]);
fpath = fullfile(files(ind).folder, files(ind).name);
data = load(fpath);
sampleScene = data.projection;
contrastScene = (sampleScene - min(sampleScene(:))) / (max(sampleScene(:)) - min(sampleScene(:)));

% Get the size of the scene database
[numY,numX] = size(sampleScene);

% Get the spatial resolution of the scenes in degrees per pixel
sceneResX = 360 / numX;


gammaScene = imadjust(contrastScene, [], [0,1], 0.5);
x = (0:numX-1) * sceneResX;
y = (-floor(numY/2):floor(numY/2)) * sceneResX;


v = 45;

t = (0:1e-3:2)';
indY = randi([1 numY]);
indX = randi([2 numX-1]);

% Convert photoreceptor spacing to pixels
prSpacingPx = 5 / sceneResX;

% Compute the offsets at each timestep (in pixels)
dx = -t * v / sceneResX + indX;

% Compute the indices for each input with circular boundary
% conditions
idx1 = 1 + mod(round(dx - prSpacingPx), numX);
idx2 = 1 + mod(round(dx), numX);
idx3 = 1 + mod(round(dx + prSpacingPx), numX);

prInputs = [gammaScene(indY,idx1)',gammaScene(indY,idx2)',gammaScene(indY,idx3)'];


MakeFigure;
imagesc(x,y,gammaScene);
hold on;
plot(x(idx1(1)), y(indY), 'o', 'markersize', 20,'linewidth',2);
plot(x(idx2(1)), y(indY), 'o', 'markersize', 20,'linewidth',2);
plot(x(idx3(1)), y(indY), 'o', 'markersize', 20,'linewidth',2);

axis('equal','tight');
cbar = colorbar;
caxis([0 1]);
colormap(gray(2^16));
xlim([0 360]);


MakeFigure;
imagesc(t, -1:1, prInputs');
cbar = colorbar;
caxis([0 1]);
xlabel('time (s)');
ylabel('photoreceptor location');
yticks(-1:1);
colormap(gray(2^16));
