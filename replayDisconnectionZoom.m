function replayDisconnectionZoom(folder,framePause)
% Replay exact saved density crops with fixed axes and initial colour limits.
if nargin < 1 || isempty(folder)
    folder = uigetdir(pwd,'Select a simulation folder or its density_zoom folder');
    if isequal(folder,0), return; end
end
if nargin < 2, framePause = 0.1; end
if isfolder(fullfile(folder,'density_zoom'))
    folder = fullfile(folder,'density_zoom');
end
files = dir(fullfile(folder,'density_step*.mat'));
assert(~isempty(files),'No density crops found in %s.',folder);
[~,order] = sort({files.name}); files = files(order);
first = load(fullfile(folder,files(1).name),'roi'); r=first.roi;
fig = figure('Name','Disconnection density replay','Color','w', ...
    'Position',[100 100 1400 550]);
ax = axes('Parent',fig);
im = imagesc(ax,r.x/r.latticeUnit,(r.y-r.flatBoundaryY)/r.latticeUnit,r.density.');
set(ax,'YDir','normal'); axis(ax,'image');
colormap(ax,gray(256)); colorbar(ax); caxis(ax,r.densityLimits);
xlabel(ax,'X/a'); ylabel(ax,'(Y-Y_{flat})/a');
for j=1:numel(files)
    if ~isgraphics(fig), break; end
    data = load(fullfile(folder,files(j).name),'roi'); r=data.roi;
    set(im,'CData',r.density.');
    title(ax,sprintf('Raw density | step %d | t=%.4g | %s',r.step,r.time,r.phase));
    drawnow;
    if framePause>0, pause(framePause); end
end
end
