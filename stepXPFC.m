function [ fields, status ] = stepXPFC( params, fields, operators, status )
% stepXPFC calculates the non-linear operators for the evolution of n and c
% fields and the performs a single time-step.

%% Helper functions for non-linear operators
% Calculate polynomial interpolation functions for the correlation
% function and its derivative with respect to c.
X_A = 1-3*fields.c.^2+2*fields.c.^3 ;
X_A_c = -6*fields.c.*(1-fields.c);

% Free energy of mixing and its derivative
F_mix = params.omega*(fields.c.*log(fields.c/params.c0)+(1-fields.c).*log( (1-fields.c)/(1-params.c0)) );
F_mix_c = params.omega*(log(fields.c/params.c0)-log((fields.c-1)/(params.c0-1) ));
% Calculation of excess terms
C_2AAn = ifftn(operators.C_2AAk.*fields.n_k);
C_2BBn = ifftn(operators.C_2BBk.*fields.n_k);
C_2n = X_A.*C_2AAn + (1-X_A).*C_2BBn;
C_2_cn = X_A_c.*C_2AAn - X_A_c.*C_2BBn;

%% Noise
% Zero forcing when either noise flag is disabled.
noise_n_k = zeros(size(fields.n_k), "like", fields.n_k);
noise_c_k = zeros(size(fields.c_k), "like", fields.c_k);
if params.noiseApp_n == true
    X = params.NA_n*randn(params.Nx,params.Ny,"gpuArray");
    Y = params.NA_n*randn(params.Nx,params.Ny,"gpuArray");

    noise_n = (Y - circshift(Y,[-1,0]) + X - circshift(X,[0,-1]))/params.dx; % This scheme assumes dx=dy
    noise_n_k = operators.filter_n.*fftn(noise_n);
end
if params.noiseApp_c == true
    X = params.NA_c*randn(params.Nx,params.Ny,"gpuArray");
    Y = params.NA_c*randn(params.Nx,params.Ny,"gpuArray");
    noise_c = (Y - circshift(Y,[-1,0]) + X - circshift(X,[0,-1]))/params.dx;
    noise_c_k = operators.filter_c.*fftn(noise_c);
end
%% Non-linear terms
NL_n = params.M_n*operators.k2.*(fftn(params.eta*fields.n.^2/2-params.chi*fields.n.^3/3 - F_mix + C_2n  ) ) + noise_n_k;
NL_c = -params.M_c*operators.k2.*fftn(-0.5*fields.n.*C_2_cn + (1 + fields.n).*F_mix_c) + noise_c_k;

%% Update step
fields.c_k = (fields.c_k+params.dt*NL_c)./(1-params.dt*operators.L_c);
% Backward second difference for n_tt and backward difference for damping.
% n_k_prev must be preserved with n_k when saving a restart.
n_k_old = fields.n_k;
fields.n_k = ((2 + params.beta_n*params.dt)*n_k_old ...
    - fields.n_k_prev + params.dt^2*NL_n) ...
    ./ (1 + params.beta_n*params.dt - params.dt^2*operators.L_n);
fields.n_k_prev = n_k_old;
% Switch n and c back to realspace
fields.n = real(ifftn(fields.n_k));
fields.c = real(ifftn(fields.c_k));
end
