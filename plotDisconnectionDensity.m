function state = plotDisconnectionDensity(fields,p,step,state,saveData)
% Fixed raw-density window encompassing both initial cores and the intervening GB.
% The view never depends on inferred boundary/core positions.
if nargin < 4, state = []; end
if nargin < 5, saveData = false; end
a = 2*pi/p.k_ini;
y0 = -p.Ny*p.dy/4;
if isempty(state)
    paddingX = 4*a; paddingY = 4*a;
    if isfield(p,'densityZoom')
        paddingX = p.densityZoom.horizontalPadding;
        paddingY = p.densityZoom.verticalPadding;
    end
    x = ((0:p.Nx-1)'-floor(p.Nx/2))*p.dx;
    y = ((0:p.Ny-1)-floor(p.Ny/2))*p.dy;
    halfWidth = p.disconnection.separation/2+paddingX;
    state.ix = find(abs(x)<=halfWidth);
    state.iy = find(y>=y0-paddingY & y<=y0+p.disconnection.stepHeight+paddingY);
    assert(~isempty(state.ix) && ~isempty(state.iy), 'Density zoom is empty.');
    state.x = x(state.ix);
    state.y = y(state.iy);
    initialDensity = gather(fields.n(state.ix,state.iy));
    state.limits = [min(initialDensity(:)) max(initialDensity(:))];
    if state.limits(2)<=state.limits(1)
        state.limits = state.limits + [-0.5 0.5];
    end
    state.figure = [];
    state.previousMeasurement = struct();
    state.trackingRows = zeros(0,10);
end
density = gather(fields.n(state.ix,state.iy));
if isempty(state.figure) || ~isgraphics(state.figure)
    state.figure = figure('Name','Disconnection pair: raw density', ...
        'Color','w','Position',[100 100 1400 550]);
    state.axes = axes('Parent',state.figure);
    state.image = imagesc(state.axes,state.x/a,(state.y-y0)/a,density.');
    set(state.axes,'YDir','normal');
    axis(state.axes,'image');
    xlim(state.axes,[state.x(1) state.x(end)]/a);
    ylim(state.axes,([state.y(1) state.y(end)]-y0)/a);
    colormap(state.axes,gray(256)); colorbar(state.axes);
    caxis(state.axes,state.limits);
    hold(state.axes,'on');
    state.boundary = plot(state.axes,NaN,NaN,'c-','LineWidth',0.7);
    state.core1 = plot(state.axes,NaN,NaN,'ro','MarkerSize',11,'LineWidth',2);
    state.core2 = plot(state.axes,NaN,NaN,'yo','MarkerSize',11,'LineWidth',2);
    xlabel(state.axes,'X/a');
    ylabel(state.axes,'(Y-Y_{flat})/a');
else
    set(state.image,'CData',density.');
end
tracking = struct('resolved',false,'x',[NaN NaN],'y',[NaN NaN], ...
    'boundaryX',[],'boundaryY',[],'measurement',struct(),'message','tracking off');
if isfield(p,'densityZoom') && isfield(p.densityZoom,'trackPositions') ...
        && p.densityZoom.trackPositions
    tracking = trackDensityZoomCores(fields,p,state.previousMeasurement);
    state.previousMeasurement = tracking.measurement;
end
set(state.boundary,'XData',tracking.boundaryX/a, ...
    'YData',(tracking.boundaryY-y0)/a);
if tracking.resolved
    set(state.core1,'XData',tracking.x(1)/a,'YData',(tracking.y(1)-y0)/a);
    set(state.core2,'XData',tracking.x(2)/a,'YData',(tracking.y(2)-y0)/a);
else
    set(state.core1,'XData',NaN,'YData',NaN);
    set(state.core2,'XData',NaN,'YData',NaN);
end
phase = 'coupled n-c';
if isfield(p,'freezeCompositionSteps') && step<=p.freezeCompositionSteps
    phase = 'composition frozen';
end
title(state.axes,sprintf('Density | step %d | t=%.4g | %s | %s', ...
    step,step*p.dt,phase,tracking.message),'Interpreter','none');
drawnow;
if saveData && isfield(p,'densityZoom') && p.densityZoom.saveData
    zoomFolder = fullfile(p.output_dir,'density_zoom');
    if ~isfolder(zoomFolder), mkdir(zoomFolder); end
    roi = struct('density',density,'x',state.x,'y',state.y, ...
        'latticeUnit',a,'flatBoundaryY',y0,'step',step,'time',step*p.dt, ...
        'densityLimits',state.limits,'phase',phase, ...
        'initialSeparation',p.disconnection.separation,'tracking',tracking);
    % Preserve every selected grid value in double precision; no image resampling.
    save(fullfile(zoomFolder,sprintf('density_step%07d.mat',step)),'roi','-v7');
    separation=NaN; height=NaN; fitRMSE=NaN;
    if isfield(tracking.measurement,'separation')
        separation=tracking.measurement.separation;
        height=tracking.measurement.height;
        fitRMSE=tracking.measurement.rmse;
    end
    state.trackingRows(end+1,:) = [step step*p.dt double(tracking.resolved) ...
        tracking.x(1) tracking.y(1) tracking.x(2) tracking.y(2) ...
        separation height fitRMSE];
    coreTable = array2table(state.trackingRows,'VariableNames', ...
        {'step','time','resolved','x1','y1','x2','y2','separation','height','fitRMSE'});
    writetable(coreTable,fullfile(zoomFolder,'core_tracking.csv'));
end
end
