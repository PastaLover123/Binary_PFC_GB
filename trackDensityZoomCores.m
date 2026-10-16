function t = trackDensityZoomCores(fields,p,previous)
% Provisional geometric core positions from the lower-GB contour.
% Failure never suppresses the raw density display.
if nargin < 3, previous = struct(); end
t = struct('resolved',false,'x',[NaN NaN],'y',[NaN NaN], ...
    'boundaryX',[],'boundaryY',[],'measurement',struct(),'message','unresolved');
try
    [~,~,~,localGB] = measureSigma5(fields,p);
    t.boundaryX = ((0:p.Nx-1)'-floor(p.Nx/2))*p.dx;
    t.boundaryY = localGB(:,1);
    m = measureDisconnectionPair(localGB,p,previous);
    t.measurement = m;
    if ~m.resolved, return; end
    t.x = [m.core1 m.core2];
    % Interpolate periodically to put the markers directly on the measured GB.
    Lx = p.Nx*p.dx; Ly = p.Ny*p.dy;
    y = unwrap(t.boundaryY*2*pi/Ly)*Ly/(2*pi);
    xx = [t.boundaryX-Lx;t.boundaryX;t.boundaryX+Lx];
    yy = [y;y;y];
    t.y = interp1(xx,yy,t.x,'linear');
    t.y = mod(t.y+Ly/2,Ly)-Ly/2;
    t.resolved = all(isfinite(t.y));
    if t.resolved, t.message = 'provisional core positions'; end
catch err
    t.resolved=false; t.x=[NaN NaN]; t.y=[NaN NaN];
    t.message = ['tracking unavailable: ' err.message];
end
end
