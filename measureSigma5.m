function [positions, profile, indicator, localGB] = measureSigma5(fields,params)
% Orientation score on the material grid (array layout [Nx,Ny]).
% Affine shear changes physical q and k together, so windows stay at reference q.
kx = gpuArray(ifftshift((-floor(params.Nx/2):ceil(params.Nx/2)-1) ...
    *2*pi/(params.Nx*params.dx))).';
ky = gpuArray(ifftshift((-floor(params.Ny/2):ceil(params.Ny/2)-1) ...
    *2*pi/(params.Ny*params.dy)));
scores = cell(1,2);
angles = [params.theta_1 params.theta_2];
width = 0.12*params.k_ini;
for grain=1:2
    t = angles(grain);
    q = params.k_ini*[cosd(t) -sind(t); sind(t) cosd(t)];
    score = zeros(size(fields.n), 'like', fields.n);
    for j=1:2
        window = exp(-((kx-q(1,j)).^2+(ky-q(2,j)).^2)/(2*width^2));
        score = score + abs(ifftn(window.*fields.n_k)).^2;
    end
    scores{grain} = score;
end
total = scores{1}+scores{2};
scale = gather(max(total(:)));
if ~isfinite(scale) || scale <= eps
    error('No crystalline signal for boundary tracking.');
end
indicator = gather((scores{1}-scores{2})./(total+1e-12*scale));
profile = mean(indicator,1).';
localGB = NaN(params.Nx,2);
if nargout >= 4
    for ix = 1:params.Nx
        localGB(ix,:) = crossings(indicator(ix,:).',params);
    end
end
positions = crossings(profile,params);
end

function positions = crossings(profile,params)
positions = [NaN NaN]; % rising boundary, falling boundary, in material Y
nextProfile = circshift(profile,-1);
rising = find(profile <= 0 & nextProfile > 0);
falling = find(profile > 0 & nextProfile <= 0);
% Do not silently pick two crossings if the orientation structure breaks up.
if numel(rising) ~= 1 || numel(falling) ~= 1, return; end
Ly = params.Ny*params.dy;
y = ((0:params.Ny-1)'-floor(params.Ny/2))*params.dy;
for b=1:2
    if b==1, j=rising; else, j=falling; end
    f = -profile(j)/(nextProfile(j)-profile(j));
    positions(b) = mod(y(j)+f*params.dy+Ly/2,Ly)-Ly/2;
end
end
