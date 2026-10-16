function [fields] = genInitialConditions(fields,params)
% genInitialConditions generates inditial condition based on a flag and
% parameters provided in params.

% Uniform values for the n and c fields
if params.initialCondition == "uniform"
    fields.n = fields.n+params.n_init;
    fields.c = fields.c+params.c_init;
% Initialize n with a mode expansion by initializing with amplitude A_ini,
% wave-vector k_ini and angles theta and a liquid barrier of
% liquidPaddingWidth around at the edge of the box.
elseif params.initialCondition == "mode expansion"
   % Generate mode expansion
   n_exp = params.n_init + rotModeExpHex2D(params.A_ini,params.k_ini,params.theta,params.Nx,params.Ny,params.dx,params.dy);
   % Mask the mode expansion and insert into the density field
   fields.n(params.liquidPaddingWidth:(params.Nx-params.liquidPaddingWidth),params.liquidPaddingWidth:(params.Ny-params.liquidPaddingWidth)) = n_exp(params.liquidPaddingWidth:(params.Nx-params.liquidPaddingWidth),params.liquidPaddingWidth:(params.Ny-params.liquidPaddingWidth));
   % smooth edges
   fields.n = imgaussfilt(real(fields.n),2,"Padding","circular");
   fields.c = fields.c + params.c_init;

% Initialize n with a bicrystal by initializing two mode expansions with
% amplitudes A_ini, wave-vectors k_ini and angles theta_1 and theta_2 with
% a liquid barrier of liquidPaddingWidth seperating them.
elseif params.initialCondition == "bicrystal"
   % Generate mode expansions
   n_1 = params.n_init + rotModeExpHex2D(params.A_ini,params.k_ini,params.theta_1,params.Nx,params.Ny,params.dx,params.dy);
   n_2 = params.n_init + rotModeExpHex2D(params.A_ini,params.k_ini,params.theta_2,params.Nx,params.Ny,params.dx,params.dy);
   % Mask the mode expansions and insert into the density field
   fields.n(params.liquidPaddingWidth:(round(params.Nx/2)-params.liquidPaddingWidth),params.liquidPaddingWidth:(params.Ny-params.liquidPaddingWidth)) = n_1(params.liquidPaddingWidth:(round(params.Nx/2)-params.liquidPaddingWidth),params.liquidPaddingWidth:(params.Ny-params.liquidPaddingWidth));
   fields.n((round(params.Nx/2)+params.liquidPaddingWidth):(params.Nx-params.liquidPaddingWidth),params.liquidPaddingWidth:(params.Ny-params.liquidPaddingWidth)) = n_2((round(params.Nx/2)+params.liquidPaddingWidth):(params.Nx-params.liquidPaddingWidth),params.liquidPaddingWidth:(params.Ny-params.liquidPaddingWidth));
   fields.c = fields.c+params.c_init;

% Initialize n with a bicrystal by initializing two mode expansions with
% different amplitudes A_iniA, A_iniB, orientations theta_1, theta_2 and
% wave-vectors chosen at the peaks of the correlation functions for pure
% A and B. The c field is then initialized such that the A and B rich
% regions have concentraiton c_A and c_B. A liquid barrier of
% liquidPaddingWidth seperates the mode expansions.
   elseif params.initialCondition == "bicrystal"
   % Generate mode expansions
   n_1 = params.n_init + rotModeExpHex2D(params.A_iniA,params.k_1A,params.theta_1,params.Nx,params.Ny,params.dx,params.dy);
   n_2 = params.n_init + rotModeExpHex2D(params.A_iniB,params.k_1B,params.theta_2,params.Nx,params.Ny,params.dx,params.dy);
   % Mask the mode expansions and insert into the density field
   fields.n(params.liquidPaddingWidth:(round(params.Nx/2)-params.liquidPaddingWidth),params.liquidPaddingWidth:(params.Ny-params.liquidPaddingWidth)) = n_1(params.liquidPaddingWidth:(round(params.Nx/2)-params.liquidPaddingWidth),params.liquidPaddingWidth:(params.Ny-params.liquidPaddingWidth));
   fields.n((round(params.Nx/2)+params.liquidPaddingWidth):(params.Nx-params.liquidPaddingWidth),params.liquidPaddingWidth:(params.Ny-params.liquidPaddingWidth)) = n_2((round(params.Nx/2)+params.liquidPaddingWidth):(params.Nx-params.liquidPaddingWidth),params.liquidPaddingWidth:(params.Ny-params.liquidPaddingWidth));
   fields.c = fields.c+params.c_init;
   % Mask the concentration and insert into the c field. Note the boundary
   % regions will be defined with c_init.
   fields.c(params.liquidPaddingWidth:(round(params.Nx/2)-params.liquidPaddingWidth),params.liquidPaddingWidth:(params.Ny-params.liquidPaddingWidth)) = params.c_A;
   fields.c((round(params.Nx/2)+params.liquidPaddingWidth):(params.Nx-params.liquidPaddingWidth),params.liquidPaddingWidth:(params.Ny-params.liquidPaddingWidth)) = params.c_B;
   fields.n = imgaussfilt(real(fields.n),2,"Padding","circular");
   fields.c = imgaussfilt(real(fields.c),10,"Padding","circular");

end
