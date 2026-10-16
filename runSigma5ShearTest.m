function runSigma5ShearTest(direction)
% direction = +1 (default), -1 (reverse shear), or 0 (unstrained control).
if nargin < 1, direction = 1; end
validateattributes(direction, {'numeric'}, {'scalar','real','finite'});
assert(ismember(direction,[-1 0 1]), 'Use direction -1, 0, or +1.');
params = genParams();
params.setNum = 3 + direction;
params.initialCondition = "sigma5";
% Two Gaussian shells support the square seed; no polynomial kernel.
params.correlationA = "square";
params.correlationB = "square";
params.k_2A = sqrt(2)*params.k_1A;
params.k_2B = sqrt(2)*params.k_1B;
params.sigma_2A = sqrt(2)*params.sigma_1A;
params.sigma_2B = sqrt(2)*params.sigma_1B;
params.c_init = 0.1;
params.n_init = -0.05;
params.k_ini = params.k_1A + params.c_init*(params.k_1B-params.k_1A);
params.theta_1 = atand(1/3);
params.theta_2 = -params.theta_1;
% Approximate minimum within the two-shell square-mode ansatz.
params.squareA10 = 0.29744;
params.squareA11 = 0.18493;
% Exactly commensurate square lattice in both material directions.
params.sigmaPeriods = 16;
L = params.sigmaPeriods*2*pi*sqrt(10)/params.k_ini;
params.Nx = 405; params.Ny = 405; % odd: no sheared Nyquist planes
params.dx = L/params.Nx; params.dy = L/params.Ny;
params.griddim = [params.Nx params.Ny];
params.gbWidth = 0.5*(2*pi/params.k_ini);
params.dt = 0.01;
params.beta_n = 1;
params.M_n = 1;
params.M_c = 0.01;
params.noiseApp_n = false; params.noiseApp_c = false;
params.gamma = 0;
params.shear.enabled = true;
params.shear.relaxSteps = 5000; % time 50 before loading
params.shear.every = 100;
params.shear.rate = direction*5e-4;
params.shear.maxStrain = 0.1;
params.maxStep = 25000; % time 250; max |gamma|=0.1
params.outStep = 1000; params.printStep = 1000;
params.trackGB = true;
params.dispOutput = true; params.makeMovie = false;
params.saveFields = true;
[~, tag] = fileparts(tempname);
params.output_dir = fullfile(pwd, sprintf('sigma5_shear_%+d_%s',direction,tag));
fprintf('Sigma5: +/- %.6f deg, L=%.6g, dx=%.6g\n', ...
    params.theta_1,L,params.dx);
fprintf('Output directory: %s\n',params.output_dir);
simulateXPFC(params);
end
