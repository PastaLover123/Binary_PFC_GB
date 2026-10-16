function p = makeDisconnectionAnnihilationParams(cMean,soluteMode)
% Small zero-shear pair-annihilation test.
if nargin < 1, cMean = 0.3; end
if nargin < 2, soluteMode = "mobile"; end
soluteMode = string(soluteMode);
assert(any(soluteMode==["mobile","frozen"]), 'Use mobile or frozen.');
p = makeSigma5DisconnectionParams(0,cMean,true);
a = 2*pi/p.k_ini;
p.Nx = 512; p.Ny = 512;
p.sigmaPeriods = 20;
L = p.sigmaPeriods*sqrt(10)*a;
p.dx = L/p.Nx; p.dy = L/p.Ny;
p.griddim = [p.Nx p.Ny];
p.disconnection.separation = 1.5*L/4;
p.shear.enabled = false;
p.shear.rate = 0;
p.maxStep = 10000;
p.preparationSteps = 100; % short density-only preparation, included in diagnostics
p.freezeCompositionSteps = p.preparationSteps;
if soluteMode=="frozen", p.freezeCompositionSteps = p.maxStep; end
p.soluteMode = soluteMode;
p.outStep = 50; p.printStep = 500;
p.fieldSaveStep = 1000;
% Raw-density inspection is independent of the provisional position tracker.
p.trackGB = false;
p.trackPair = false;
p.dispOutput = false;
p.densityZoom.enabled = true;
p.densityZoom.saveData = true;
p.densityZoom.horizontalPadding = 4*a;
p.densityZoom.verticalPadding = 4*a;
p.pairTracking.fitWidth = 1.2*a;
p.pairTracking.minHeight = 0.25*p.disconnection.stepHeight;
[~, tag] = fileparts(tempname);
p.output_dir = fullfile(pwd,sprintf('annihilation_c%.2f_%s_%s',cMean,soluteMode,tag));
end
