function m = measureSigma5Profile(fields,p)
% Gather only the Y orientation profile, not a 2D indicator or local contours.
kx=gpuArray(ifftshift((-floor(p.Nx/2):ceil(p.Nx/2)-1) ...
    *2*pi/(p.Nx*p.dx))).';
ky=gpuArray(ifftshift((-floor(p.Ny/2):ceil(p.Ny/2)-1) ...
    *2*pi/(p.Ny*p.dy)));
angles=[p.theta_1 p.theta_2]; scores=cell(1,2);
for grain=1:2
    theta=angles(grain);
    Q=p.k_ini*[cosd(theta) -sind(theta);sind(theta) cosd(theta)];
    score=zeros(p.griddim,'like',fields.n);
    for j=1:2
        W=exp(-((kx-Q(1,j)).^2+(ky-Q(2,j)).^2)/(2*(0.12*p.k_ini)^2));
        score=score+abs(ifftn(fields.n_k.*W)).^2;
    end
    scores{grain}=score;
end
total=scores{1}+scores{2};
scale=max(total(:));
profile=gather(mean((scores{1}-scores{2})./(total+max(1e-12*scale,eps)),1)).';
m=struct('positions',[NaN NaN],'profile',profile, ...
    'strength',mean(abs(profile)));
if any(~isfinite(profile)) || m.strength<0.25, return; end
next=circshift(profile,-1);
rising=find(profile<=0 & next>0);
falling=find(profile>0 & next<=0);
if numel(rising)~=1 || numel(falling)~=1, return; end
Ly=p.Ny*p.dy;
y=((0:p.Ny-1)'-floor(p.Ny/2))*p.dy;
for b=1:2
    if b==1, j=rising; else, j=falling; end
    f=-profile(j)/(next(j)-profile(j));
    m.positions(b)=mod(y(j)+f*p.dy+Ly/2,Ly)-Ly/2;
end
end
