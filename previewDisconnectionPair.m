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
[~,~,~,localGB] = measureSigma5(fields,p);
a = 2*pi/p.k_ini;
x = ((0:p.Nx-1)-floor(p.Nx/2))*p.dx/a;
y = ((0:p.Ny-1)-floor(p.Ny/2))*p.dy/a;
n = gather(fields.n);
yc = -p.Ny*p.dy/(4*a);
cores = [-1 1]*p.disconnection.separation/(2*a);
figure('Name','Disconnection initializer preview: no evolution','Color','w');
tiledlayout(2,2);
nexttile;
imagesc(x,y,n.'); axis image; set(gca,'YDir','normal');
xlim([cores(1)-4 cores(2)+4]); ylim(yc+[-4 5]);
hold on; plot(x,localGB(:,1)/a,'r-','LineWidth',1.2);
title('Seeded terrace and density'); xlabel('X/a'); ylabel('Y/a');
for j=1:2
    nexttile;
    imagesc(x,y,n.'); axis image; set(gca,'YDir','normal');
    xlim(cores(j)+[-3 3]); ylim(yc+[-3 4]);
    hold on; plot(x,localGB(:,1)/a,'r-','LineWidth',1.2);
    title(sprintf('Core %d near X/a=%.2f',j,cores(j)));
    xlabel('X/a'); ylabel('Y/a');
end
colormap(gray);
nexttile;
plot(x,localGB(:,1)/a-yc,'k-','LineWidth',1.2); hold on;
yline(0,'--'); yline(p.disconnection.stepHeight/a,'--');
xline(cores(1),':'); xline(cores(2),':');
xlabel('X/a'); ylabel('(Y_{GB}-Y_{flat})/a');
title('Local lower-GB contour');
fprintf('Preview only: no timesteps performed.\n');
fprintf('Prescribed h/a=%.5f, b/a=%.5f, core separation/a=%.5f\n', ...
    p.disconnection.stepHeight/a,p.disconnection.burgers/a, ...
    p.disconnection.separation/a);
fprintf('Run this exact configuration with runDisconnectionAnnihilationTest(p).\n');
end
