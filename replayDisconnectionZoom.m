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
xlim(ax,[r.x(1) r.x(end)]/r.latticeUnit);
ylim(ax,([r.y(1) r.y(end)]-r.flatBoundaryY)/r.latticeUnit);
colormap(ax,gray(256)); colorbar(ax); caxis(ax,r.densityLimits);
hold(ax,'on');
boundary = plot(ax,NaN,NaN,'c-','LineWidth',0.7);
core1 = plot(ax,NaN,NaN,'ro','MarkerSize',11,'LineWidth',2);
core2 = plot(ax,NaN,NaN,'yo','MarkerSize',11,'LineWidth',2);
xlabel(ax,'X/a'); ylabel(ax,'(Y-Y_{flat})/a');
for j=1:numel(files)
    if ~isgraphics(fig), break; end
    data = load(fullfile(folder,files(j).name),'roi'); r=data.roi;
    set(im,'CData',r.density.');
    set(boundary,'XData',NaN,'YData',NaN);
    set(core1,'XData',NaN,'YData',NaN);
    set(core2,'XData',NaN,'YData',NaN);
    if isfield(r,'tracking')
        t=r.tracking; a=r.latticeUnit; y0=r.flatBoundaryY;
        set(boundary,'XData',t.boundaryX/a,'YData',(t.boundaryY-y0)/a);
        if t.resolved
            set(core1,'XData',t.x(1)/a,'YData',(t.y(1)-y0)/a);
            set(core2,'XData',t.x(2)/a,'YData',(t.y(2)-y0)/a);
        end
    end
    title(ax,sprintf('Raw density | step %d | t=%.4g | %s',r.step,r.time,r.phase));
    drawnow;
    if framePause>0, pause(framePause); end
end
end
