
% Creates a mode expansion for a hexagonal phase given an amplitude, input
% k-vector and an angle using gpu arrays.
function[n] = rotModeExpHex2D(A,k_in,theta,Nx,Ny,dx,dy)
% rotModeExpHex2D generates a two-dimensional mode-expansion for a
% hexagonal phase with amplitude A, wave-vector k_in and angle theta.
x = dx*gpuArray.linspace(1,Nx,Nx);
y = dy*gpuArray.linspace(1,Ny,Ny);

[Y0, X0] = meshgrid(y,x);
XY = [X0(:) Y0(:)];
R=[cosd(theta) -sind(theta); sind(theta) cosd(theta)];
rotXY=XY*R';
X = reshape(rotXY(:,1), size(X0,1), []);
Y = reshape(rotXY(:,2), size(Y0,1), []);

n = 2*A*cos(k_in*X) + 4*A*cos(k_in*X/2).*cos(sqrt(3)/2*k_in*Y);
