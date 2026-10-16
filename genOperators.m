function [operators] = genOperators(params)
% genOperators defines and returns the operators structure which contains the
% density-density 2 point correlation functions in Fourier space for pure A,
% and pure B species. Other relevant Fourier space operators that do not
% shift dynamically are also contained within this structure
operators = struct();
% Calculate k^2 and k^4 operator for the linear components of n and c. Only
% k2 is required for the non-linear terms.
dkx = 2*pi/(params.Nx*params.dx);
dky = 2*pi/(params.Ny*params.dy);
kx = ifftshift(dkx*gpuArray.linspace(-params.Nx/2,params.Nx/2-1,params.Nx));
ky = ifftshift(dky*gpuArray.linspace(-params.Ny/2,params.Ny/2-1,params.Ny));
[kY,kX] = meshgrid(ky,kx);
operators.k2 = kX.^2+kY.^2;
% Locally calculate the k operator for the correlation functions as this is
% not used in the dynamics.
k = sqrt(operators.k2);
%operators.k4 = k2.^2;

% Define the Correlation functions
if params.correlationA == "FCC" || params.correlationA == "DC"
    operators.C_2AAk = max(exp(-params.T/params.T0_1A).*exp( -(k-params.k_1A).^2/2/params.sigma_1A^2 ),exp(-params.T/params.T0_2A).*exp( -(k-params.k_2A).^2/2/params.sigma_2A^2 ));
elseif params.correlationA == "BCC"
    operators.C_2AAk = exp(-params.T/params.T0_1A).*exp( -(k-params.k_1A).^2/2/params.sigma_1A^2 );
end

if params.correlationB == "FCC" || params.correlationB == "DC"
    operators.C_2BBk = max(exp(-params.T/params.T0_1B).*exp( -(k-params.k_1B).^2/2/params.sigma_1B^2 ),exp(-params.T/params.T0_2B).*exp( -(k-params.k_2B).^2/2/params.sigma_2B^2 ));
elseif params.correlationB == "BCC"
    operators.C_2BBk = exp(-params.T/params.T0_1B).*exp( -(k-params.k_1B).^2/2/params.sigma_1B^2 );
end

% Define linear operators for the n and c fields
operators.L_n = -params.M_n*operators.k2;
operators.L_c = -params.M_c*params.W_c*operators.k2.^2;
% Wavelength cutoff filter for the noise to remove fluctuations on scales
% smaller than the lattice constant.
if params.correlationB == "FCC" || params.correlationB == "DC"
    operators.filter_n = (sqrt(operators.k2) < params.k_2B);
    operators.filter_c = (sqrt(operators.k2) < params.k_2B);
elseif params.correlationB == "BCC"
    operators.filter_n = (sqrt(operators.k2) < params.k_1B);
    operators.filter_c = (sqrt(operators.k2) < params.k_1B);
end
