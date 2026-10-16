function m = measureDisconnectionPair(localGB,p,previous)
% Fit a periodic terrace to the lower-GB contour. This is a geometric proxy;
% loss of a resolved terrace is not proof of topological annihilation.
if nargin < 3, previous = struct(); end
m = struct('centre',NaN,'fittedSeparation',NaN,'separation',NaN, ...
    'height',NaN,'rmse',NaN,'resolved',false,'core1',NaN,'core2',NaN);
y = localGB(:,1);
if any(~isfinite(y)), return; end
Ly = p.Ny*p.dy;
y = unwrap(y*2*pi/Ly)*Ly/(2*pi);
L = p.Nx*p.dx;
x = ((0:p.Nx-1)'-floor(p.Nx/2))*p.dx;
w = p.pairTracking.fitWidth;
guess = [0 p.disconnection.separation];
if isfield(previous,'centre') && isfinite(previous.centre) ...
        && isfinite(previous.fittedSeparation)
    guess = [previous.centre previous.fittedSeparation];
end
opt = optimset('Display','off','MaxIter',150,'TolX',1e-6,'TolFun',1e-10);
fit = fminsearch(@cost,guess,opt);
if fit(2) <= w || fit(2) >= L/2, return; end
[~,ab,residual] = cost(fit);
m.centre = mod(fit(1)+L/2,L)-L/2;
m.fittedSeparation = fit(2);
m.height = ab(2); m.rmse = sqrt(mean(residual.^2));
m.resolved = m.height >= p.pairTracking.minHeight ...
    && fit(2) > 2*w && m.rmse < 0.35*m.height;
if m.resolved
    m.separation = fit(2);
    m.core1 = mod(m.centre-fit(2)/2+L/2,L)-L/2;
    m.core2 = mod(m.centre+fit(2)/2+L/2,L)-L/2;
end

    function [value,ab,residual] = cost(z)
        if z(2) <= w || z(2) >= L/2
            value = 1e6+sum(z.^2); ab=[NaN;NaN]; residual=NaN;
            return;
        end
        angle = pi*z(2)/L;
        window = 0.5*(1+tanh((cos(2*pi*(x-z(1))/L)-cos(angle)) ...
            /(sin(angle)*2*pi*w/L)));
        design = [ones(size(x)) window];
        ab = design\y;
        residual = y-design*ab;
        value = mean(residual.^2);
    end
end
