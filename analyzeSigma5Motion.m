function T = analyzeSigma5Motion(raw,p)
% raw columns: step,time,gamma,y1,y2,meanN,meanC,orientationStrength.
% Once an identity is lost, do not invent a periodic trajectory across the gap.
n=size(raw,1); Ly=p.Ny*p.dy;
unwrapped=NaN(n,2);
for b=1:2
    if ~isfinite(raw(1,3+b)), continue; end
    unwrapped(1,b)=raw(1,3+b);
    for j=2:n
        if ~isfinite(raw(j,3+b)) || ~isfinite(unwrapped(j-1,b)), break; end
        delta=mod(raw(j,3+b)-raw(j-1,3+b)+Ly/2,Ly)-Ly/2;
        unwrapped(j,b)=unwrapped(j-1,b)+delta;
    end
end
timeSinceShear=raw(:,2)-p.shear.relaxSteps*p.dt;
baseline=find(raw(:,1)==p.shear.relaxSteps,1);
displacement=NaN(n,2);
if ~isempty(baseline)
    displacement=unwrapped-unwrapped(baseline,:);
end
velocity=NaN(n,2);
for j=1:n
    if timeSinceShear(j)<0, continue; end
    ii=max(1,j-p.velocityWindow+1):j;
    ii=ii(timeSinceShear(ii)>=0);
    if numel(ii)<3, continue; end
    for b=1:2
        if all(isfinite(displacement(ii,b)))
            fit=polyfit(timeSinceShear(ii),displacement(ii,b),1);
            velocity(j,b)=fit(1);
        end
    end
end
T=array2table([raw(:,1:5) timeSinceShear displacement velocity abs(velocity) ...
    raw(:,6:8)],'VariableNames', ...
    {'step','time','gamma','y1','y2','timeSinceShear', ...
    'displacement1','displacement2','velocity1','velocity2','speed1','speed2', ...
    'meanN','meanC','orientationStrength'});
end
