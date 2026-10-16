function runDisconnectionAnnihilationTest(input,soluteMode)
% Pass the structure returned by previewDisconnectionPair to run that exact seed.
if nargin < 1, input = 0.3; end
if nargin < 2, soluteMode = "mobile"; end
if isstruct(input), p = input;
else, p = makeDisconnectionAnnihilationParams(input,soluteMode);
end
fprintf('Zero-shear pair annihilation: %dx%d, c=%.3f, solute=%s\n', ...
    p.Nx,p.Ny,p.c_init,p.soluteMode);
fprintf('Density-only preparation: %d steps; total: %d steps.\n', ...
    p.preparationSteps,p.maxStep);
fprintf('Output: %s\n',p.output_dir);
simulateXPFC(p);
end
