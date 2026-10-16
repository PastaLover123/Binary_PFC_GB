function [p,fields] = previewDisconnectionPair(input)
% Initialize and visualize only: no time integration and no output files.
% p=previewDisconnectionPair(); runDisconnectionAnnihilationTest(p);
if nargin < 1, input = 0.3; end
if isstruct(input), p = input;
else, p = makeDisconnectionAnnihilationParams(input);
end
fields.n = zeros(p.griddim,'gpuArray');
fields.c = zeros(p.griddim,'gpuArray');
fields = genInitialConditions(fields,p);
op = genOperators(p);
fields.n_k = fftn(fields.n).*op.spectralMask;
fields.c_k = fftn(fields.c).*op.spectralMask;
fields.n_k_prev = fields.n_k;
fields.n = real(ifftn(fields.n_k)); fields.c = real(ifftn(fields.c_k));
% No orientation score, fitted contour or core-position detection is used.
plotDisconnectionDensity(fields,p,0,[],false);
a = 2*pi/p.k_ini;
fprintf('Preview only: no timesteps performed.\n');
fprintf('Prescribed h/a=%.5f, b/a=%.5f, initial separation/a=%.5f\n', ...
    p.disconnection.stepHeight/a,p.disconnection.burgers/a, ...
    p.disconnection.separation/a);
fprintf('Run this exact configuration with runDisconnectionAnnihilationTest(p).\n');
end
