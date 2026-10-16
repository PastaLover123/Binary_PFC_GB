function testDisconnectionPairTracker()
% CPU-only geometric tests; this does not identify topological core content.
p = makeDisconnectionAnnihilationParams();
a = 2*pi/p.k_ini; L = p.Nx*p.dx; Ly = p.Ny*p.dy;
x = ((0:p.Nx-1)'-floor(p.Nx/2))*p.dx;
centre = 3*a; separation = p.disconnection.separation;
angle = pi*separation/L; w = p.pairTracking.fitWidth;
window = 0.5*(1+tanh((cos(2*pi*(x-centre)/L)-cos(angle)) ...
    /(sin(angle)*2*pi*w/L)));
localGB = [-Ly/4+p.disconnection.stepHeight*window, ...
    Ly/4+zeros(size(x))];
m = measureDisconnectionPair(localGB,p);
assert(m.resolved);
assert(abs(m.separation-separation)<1e-3*a);
assert(abs(m.centre-centre)<1e-3*a);
assert(abs(m.height-p.disconnection.stepHeight)<1e-3*a);
localGB(:,1) = -Ly/4;
m = measureDisconnectionPair(localGB,p);
assert(~m.resolved && isnan(m.separation));
localGB(1,1) = NaN;
m = measureDisconnectionPair(localGB,p);
assert(~m.resolved);
disp('Known terrace, translated centre, flat boundary and invalid contour passed.');
end
