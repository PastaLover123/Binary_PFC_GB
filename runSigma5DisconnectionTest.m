function runSigma5DisconnectionTest(direction, cMean, withPair)
% 1024-square test: density-only relaxation, solute relaxation, then shear.
% Examples: (), (0), (-1), (1,0.3,false) for a flat-boundary control.
if nargin < 1, direction = 1; end
if nargin < 2, cMean = 0.3; end
if nargin < 3, withPair = true; end
params = makeSigma5DisconnectionParams(direction,cMean,withPair);
L = params.Nx*params.dx; a = 2*pi/params.k_ini;
[~, tag] = fileparts(tempname);
params.output_dir = fullfile(pwd,sprintf( ...
    'sigma5_pair%d_c%.2f_shear%+d_%s',logical(withPair),cMean,direction,tag));
fprintf('1024x1024, L=%.6g, dx=%.6g (%.3f points/lattice unit)\n', ...
    L,params.dx,a/params.dx);
fprintf('Seeded pair: h=%.6g, b=%.6g, separation=%.6g; c=%.3f\n', ...
    params.disconnection.stepHeight,params.disconnection.burgers, ...
    params.disconnection.separation,cMean);
fprintf('Steps 1:5000: n only; 5001:10000: n+c; 10001:30000: shear.\n');
fprintf('Output directory: %s\n',params.output_dir);
simulateXPFC(params);
end
