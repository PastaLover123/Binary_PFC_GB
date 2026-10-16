function [params] = genParams()
% genParams defines and returns the params structure which contains all
% parameters that are not dynamically evolving over the course of the
% simulation.

params = struct();
%The dataset number, also serves to name the set and determine seed
params.setNum=0001;

%% Phase Diagram / Misc. Model Parameters
% Expansion point for the c variable
params.c0 = 0.5;
% Ideal free energy polynomial prefactors
params.eta = 1.4;
params.chi = 1;
% Mixing free energy prefactor
params.omega = 0.002;
% Interface width for c-field (nabla ^2 prefactor)
params.W_c = 0.5;

%% Correlation Function Parameters
% Correlation type
params.correlationA = "BCC";
params.correlationB = "BCC";

% Lattice constants
params.a1 = 1;
params.a2 = 0.95;
% Wave-vectors
if params.correlationA == "FCC"
    params.k_1A = 2*pi/params.a1;
    params.k_2A = sqrt(2)*2*pi/params.a1;
elseif params.correlationA == "BCC"
    params.k_1A = 2*pi/params.a1;
elseif params.correlationA == "DC"
    params.k_1A = 2*pi/params.a1;
    params.k_2A = sqrt(8)/sqrt(3)*2*pi/params.a1;
end
if params.correlationB == "FCC"
    params.k_1B = 2*pi/params.a2;
    params.k_2B = sqrt(2)*2*pi/params.a2;
elseif params.correlationB == "BCC"
    params.k_1B = 2*pi/params.a2;
elseif params.correlationB == "DC"
    params.k_1B = 2*pi/params.a2;
    params.k_2B = sqrt(8)/sqrt(3)*2*pi/params.a2;
end

% Debye-Waller prefactor (e^(-T/T0) temperature scaling) for each
% correlation peak.
params.T0_1A = 100;
params.T0_1B = 100;
params.T0_2A = 100;
params.T0_2B = 100;
% Elastic (Gaussian) width of each wave-vector peak
params.sigma_1A = 1;
params.sigma_1B = 1;

if params.correlationA == "FCC"
    params.sigma_2A = sqrt(2)*params.sigma_1A;
elseif params.correlationA == "DC"
    params.sigma_2A = params.sigma_1A;
end

if params.correlationB == "FCC"
    params.sigma_2B = sqrt(2)*params.sigma_1B;
elseif params.correlationA == "DC"
    params.sigma_2B = params.sigma_1B;
end


%% Discretization/System Size
% System size in x and y direction
params.Nx = 512;
params.Ny = 1024;
params.griddim = [params.Nx params.Ny];
% Discretization in x and y direction
params.dx = params.a2/8;
params.dy = params.dx;

%% Dynamics
% Discretization in time
params.dt = 1;
% Density: n_tt + beta_n*n_t = M_n*laplacian(mu_n) + noise_n.
% Composition: c_t = M_c*laplacian(mu_c) + noise_c.
% beta_n controls damping; long-time density mobility is M_n/beta_n.
params.beta_n = 1;
% Mobility of n and c fields
params.M_n = 1;
params.M_c = 1;
% Timestep parameters, maxStep should be a multiple of outStep
params.maxStep = 5000;
params.outStep = 100;
params.printStep = params.outStep;

%% Initial Conditions
% Flag for initial condition :
% "uniform", "mode expansion", "bicrystal", "bicrystal A B"
params.initialCondition = "bicrystal";

% Mean field initial conditions for n,c
params.c_init = 0.1;
params.n_init=0.0;

% Isothermal quench temperature
params.T = 1;

% padding width for mode expansions
params.liquidPaddingWidth = 8;

% Initial wave-vector for bicrystal
params.k_ini = params.k_1A + 0.25*(params.k_1B - params.k_1A);
% Input amplitude for mode-expansion
params.A_ini = 0.25;
% Mode expansion angle
params.theta=0;

% Input amplitudes for two-sided concentration bicrystal
params.A_iniA = 0.25;
params.A_iniB = 0.25;

% Bicrystal mode expansion with differing angles
params.theta_1 = 4;
params.theta_2 = -4;

% Concentration of A and B rich regions for two-sided concentration
% bicrystal
params.c_A = 0.1;
params.c_B = 0.9;


%% I/O
% Set this to false for remote runs to avoid cluttered output
params.dispOutput = true;
% Flag for generation of low-quality movie
params.makeMovie = true;
% Frames per second for movie
params.movieFramerate = 5;
% Set output directory
params.output_dir = pwd;

%% Noise
% Boolean variables defining whether or not noise is applied dynamicall for
% n and c fields.
params.noiseApp_n = true;
params.noiseApp_c = true;
% Density noise amplitude numerical prefactor
params.NAP_n = 0.001;
% Density noise amplitude with appropriate scaling
params.NA_n = params.NAP_n*params.M_n*params.T/sqrt(params.dx*params.dy*params.dt);
% Composition noise amplitude numerical prefactor
params.NAP_c = 0.001;
% Composition noise amplitude with appropriate scaling
params.NA_c = params.NAP_c*params.M_c*params.T/sqrt(params.dx*params.dy*params.dt);
end
