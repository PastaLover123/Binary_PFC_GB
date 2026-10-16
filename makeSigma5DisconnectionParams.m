function params = makeSigma5DisconnectionParams(direction,cMean,withPair)
% Parameters only: callers may adjust duration/output before simulateXPFC.
if nargin < 1, direction = 1; end
if nargin < 2, cMean = 0.3; end
if nargin < 3, withPair = true; end
params = makeSigma5ShearParams(direction);
validateattributes(cMean,{'numeric'},{'scalar','real','>',0,'<',1});
validateattributes(withPair,{'logical','numeric'},{'scalar'});
params.c_init = cMean;
params.M_c = 0.1; % ten times the small test's solute mobility after release
% Weighted first-peak estimate follows the actual composition interpolant.
fractionB = 3*cMean^2-2*cMean^3;
params.k_ini = params.k_1A + fractionB*(params.k_1B-params.k_1A);
a = 2*pi/params.k_ini;
params.Nx = 1024; params.Ny = 1024;
params.sigmaPeriods = 40;
L = params.sigmaPeriods*sqrt(10)*a;
params.dx = L/params.Nx; params.dy = L/params.Ny;
params.griddim = [params.Nx params.Ny];
params.gbWidth = 0.5*a;
params.freezeCompositionSteps = 5000; % n relaxes alone until time 50
params.shear.relaxSteps = 10000; % then n and c relax until time 100
params.maxStep = 30000; % shear to |gamma|=0.1 from time 100 to 300
% Square-Sigma5 crystallographic construction:
% t_plus=(-a,3a)/sqrt(10), t_minus=(a,3a)/sqrt(10) are
% lattice translations of the +/- grains.
% b=t_minus-t_plus=(2a/sqrt(10),0), h=3a/sqrt(10), b/h=2/3.
% A terrace at +h with grain 1 shifted by +b has the same local registry,
% up to common translation t_minus, as the flat boundary.
params.disconnection.enabled = logical(withPair);
params.disconnection.burgers = 2*a/sqrt(10);
params.disconnection.stepHeight = 3*a/sqrt(10);
params.disconnection.separation = L/4;
params.disconnection.coreWidth = 0.5*a;
params.disconnection.slipDecayLength = 4*a;
end
