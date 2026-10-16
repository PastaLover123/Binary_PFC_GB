function testDensityDynamics()
% CPU-only checks of the actual step routine; no GPU/toolbox needed.
params = struct('dt', 0.01, 'beta_n', 2, 'M_n', 1, 'M_c', 1, ...
    'eta', 0, 'chi', 0, 'omega', 0, 'c0', 0.5, ...
    'noiseApp_n', false, 'noiseApp_c', false);
[Nx, Ny] = deal(16, 24);
[kx, ky] = ndgrid([0:7 -8:-1], [0:11 -12:-1]);
operators.k2 = kx.^2 + ky.^2;
operators.C_2AAk = zeros(Nx, Ny);
operators.C_2BBk = zeros(Nx, Ny);
operators.L_n = -params.M_n*operators.k2;
operators.L_c = -0.5*params.M_c*operators.k2.^2;
[x, ~] = ndgrid((0:Nx-1)'/Nx, (0:Ny-1)/Ny);
fields.n = 0.1 + 0.01*cos(2*pi*x);
fields.c = 0.4 + zeros(Nx, Ny);
fields.n_k = fftn(fields.n);
fields.c_k = fftn(fields.c);
fields.n_k_prev = fields.n_k;
initial = fields;
status = struct();
[fields, ~] = stepXPFC(params, fields, operators, status);
expected = (1 + params.beta_n*params.dt)*initial.n_k ...
    ./ (1 + params.beta_n*params.dt - params.dt^2*operators.L_n);
assert(max(abs(fields.n_k(:)-expected(:))) < 1e-12);
assert(isequal(fields.n_k_prev, initial.n_k));
for j = 1:100
    [fields, ~] = stepXPFC(params, fields, operators, status);
end
assert(abs(mean(fields.n(:))-mean(initial.n(:))) < 1e-12);
assert(abs(mean(fields.c(:))-mean(initial.c(:))) < 1e-12);

% Isolate mixing: its positive chemical potential produces -M_n*k^2*F_mix.
fields = initial;
fields.n(:) = 0;
fields.n_k = fftn(fields.n);
fields.n_k_prev = fields.n_k;
fields.c = 0.4 + 0.02*cos(2*pi*x);
fields.c_k = fftn(fields.c);
params.omega = 0.002;
g = params.omega*(fields.c.*log(fields.c/params.c0) ...
    + (1-fields.c).*log((1-fields.c)/(1-params.c0)));
expected = -params.dt^2*params.M_n*operators.k2.*fftn(g) ...
    ./ (1 + params.beta_n*params.dt - params.dt^2*operators.L_n);
[fields, ~] = stepXPFC(params, fields, operators, status);
assert(max(abs(fields.n_k(:)-expected(:))) < 1e-12);
% Composition is exactly frozen, including its Fourier state, then released.
params.freezeCompositionSteps = 2;
status.timestep = 1;
params.noiseApp_c = true; % even enabled noise must be skipped while frozen
before = fields;
[fields, ~] = stepXPFC(params,fields,operators,status);
assert(isequal(fields.c,before.c) && isequal(fields.c_k,before.c_k));
status.timestep = 3;
params.noiseApp_c = false;
[released, ~] = stepXPFC(params,fields,operators,status);
assert(max(abs(released.c(:)-fields.c(:))) > 0);

% Projected even-grid shear keeps Fourier and real fields consistent.
operators.spectralMask = ones(Nx,Ny);
operators.spectralMask(Nx/2+1,:) = 0;
operators.spectralMask(:,Ny/2+1) = 0;
operators.k2 = kx.^2+(ky-0.1*kx).^2;
operators.L_n = -params.M_n*operators.k2;
operators.L_c = -0.5*params.M_c*operators.k2.^2;
fields = initial;
% Deliberately inject a Hermitian pair on a forbidden Nyquist plane.
fields.n_k(Nx/2+1,2) = 1+2i;
fields.n_k(Nx/2+1,Ny) = 1-2i;
fields.n_k_prev = fields.n_k;
fields.n = real(ifftn(fields.n_k));
[fields, ~] = stepXPFC(params,fields,operators,status);
assert(max(abs(imag(ifftn(fields.n_k))), [], 'all') < 1e-12);
assert(max(abs(fields.n_k(operators.spectralMask==0))) == 0);
assert(max(abs(fields.n_k_prev(operators.spectralMask==0))) == 0);
disp('Density, mixing sign, conservation, frozen composition and even-grid projection checks passed.');
end
