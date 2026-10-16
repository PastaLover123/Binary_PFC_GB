function p = makeDisconnectionAnnihilationParams(cMean,soluteMode,modeMultiple)
% Small zero-shear pair-annihilation test.
if nargin < 1, cMean = 0.3; end
if nargin < 2, soluteMode = "mobile"; end
if nargin < 3, modeMultiple = 2; end
validateattributes(modeMultiple,{'numeric'},{'scalar','integer','positive'});
soluteMode = string(soluteMode);
assert(any(soluteMode==["mobile","frozen"]), 'Use mobile or frozen.');
p = makeSigma5DisconnectionParams(0,cMean,true);
a = 2*pi/p.k_ini;
p.Nx = 768; p.Ny = 512;
p.sigmaPeriods = [30 20]; % rectangular CSL repeats along X and Y
Lx = p.sigmaPeriods(1)*sqrt(10)*a;
Ly = p.sigmaPeriods(2)*sqrt(10)*a;
p.dx = Lx/p.Nx; p.dy = Ly/p.Ny;
p.griddim = [p.Nx p.Ny];
% Stay below half the periodic X length; otherwise the shorter image arc wins.
p.disconnection.separation = 0.45*Lx;
p.disconnection.modeMultiple = modeMultiple;
p.disconnection.stepHeight = p.disconnection.modeMultiple*3*a/sqrt(10);
p.disconnection.burgers = p.disconnection.modeMultiple*2*a/sqrt(10);
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
p.densityZoom.trackPositions = true;
p.densityZoom.horizontalPadding = 4*a;
p.densityZoom.verticalPadding = 4*a;
p.pairTracking.fitWidth = 1.2*a;
p.pairTracking.minHeight = 0.25*p.disconnection.stepHeight;
[~, tag] = fileparts(tempname);
p.output_dir = fullfile(pwd,sprintf('annihilation_c%.2f_%s_%s',cMean,soluteMode,tag));
end
