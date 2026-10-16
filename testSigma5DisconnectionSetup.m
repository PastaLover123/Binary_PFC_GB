function testSigma5DisconnectionSetup()
% Requires MATLAB GPU; checks the full 1024-square initializer and its geometry.
p = makeSigma5DisconnectionParams();
a = 2*pi/p.k_ini;
assert(p.Nx==1024 && p.Ny==1024);
assert(abs(p.disconnection.burgers/p.disconnection.stepHeight-2/3)<1e-12);
tp = a*[-1;3]/sqrt(10); tm = a*[1;3]/sqrt(10);
for grain=1:2
    if grain==1, theta=p.theta_1; t=tp; else, theta=p.theta_2; t=tm; end
    R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
    assert(max(abs(R.'*t/a-[0;1]))<1e-12);
end
fields.n = zeros(p.griddim,'gpuArray');
fields.c = zeros(p.griddim,'gpuArray');
fields = genInitialConditions(fields,p);
fields.n_k = fftn(fields.n); fields.c_k = fftn(fields.c);
op = genOperators(p);
discard = fields.n_k.*(1-op.spectralMask);
fraction = gather(sum(abs(discard(:)).^2)/sum(abs(fields.n_k(:)).^2));
fields.n_k = fields.n_k.*op.spectralMask;
fields.n_k_prev = fields.n_k;
fields.c_k = fields.c_k.*op.spectralMask;
fields.n = real(ifftn(fields.n_k)); fields.c = real(ifftn(fields.c_k));
[pos,~,~,localGB] = measureSigma5(fields,p);
assert(all(isfinite(pos)) && all(isfinite(localGB),'all'));
assert(abs(gather(mean(fields.n(:)))-p.n_init)<1e-12);
assert(abs(gather(mean(fields.c(:)))-p.c_init)<1e-12);
x = ((0:p.Nx-1)'-floor(p.Nx/2))*p.dx;
middle = abs(x)<p.disconnection.separation/4;
outside = abs(x)>p.disconnection.separation;
height = median(localGB(middle,1))-median(localGB(outside,1));
assert(abs(height-p.disconnection.stepHeight)<0.3*a);
p.gamma = 0.1; op = genOperators(p);
before = fields.c;
status.timestep = 1;
[fields,~] = stepXPFC(p,fields,op,status);
assert(isequal(before,fields.c));
assert(gather(max(abs(imag(ifftn(fields.n_k))),[],'all'))<1e-11);
fprintf('PASS: h from contour=%.6g; prescribed h=%.6g; removed Nyquist power=%.3e\n', ...
    height,p.disconnection.stepHeight,fraction);
end
