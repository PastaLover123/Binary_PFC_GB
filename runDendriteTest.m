function runDendriteTest()
% Crystal seed in a homogeneous liquid; deterministic growth smoke test.
params = genParams();
params.setNum = 2;
params.Nx = 512;
params.Ny = 512;
params.griddim = [params.Nx params.Ny];
params.initialCondition = "crystal seed";
params.n_init = -0.05;
params.c_init = 0.1;
params.k_ini = params.k_1A + params.c_init*(params.k_1B-params.k_1A);
params.A_ini = 0.15;
params.theta = 0;
% Radius/width are physical lengths, not grid-point counts.
params.seedRadius = 5*params.a1;
params.seedWidth = params.a1;
params.dt = 0.01;
params.beta_n = 1;
params.M_n = 1;
params.M_c = 0.01;
params.maxStep = 20000;
params.outStep = 1000;
params.printStep = 1000;
params.noiseApp_n = false;
params.noiseApp_c = false;
params.dispOutput = true;
params.makeMovie = false;
params.saveFields = true;
% Unique output folder avoids overwriting earlier tests.
[~, tag] = fileparts(tempname);
params.output_dir = fullfile(pwd, ['dendrite_test_' tag]);
fprintf('Output directory: %s\n', params.output_dir);
simulateXPFC(params);
end
