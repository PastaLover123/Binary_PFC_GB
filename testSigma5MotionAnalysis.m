function testSigma5MotionAnalysis()
% CPU-only checks: opposite velocities, periodic crossing, missing samples.
p=makeSigma5CompositionParams(0.3);
p.shear.relaxSteps=200; p.velocityWindow=5;
time=(0:1:12).'; step=time/p.dt; Ly=p.Ny*p.dy;
y1=Ly/2-2+0.5*time; y2=-Ly/4-0.3*time;
wrap=@(y) mod(y+Ly/2,Ly)-Ly/2;
raw=[step time zeros(size(time)) wrap(y1) wrap(y2) ...
    -0.05+zeros(size(time)) 0.3+zeros(size(time)) ones(size(time))];
T=analyzeSigma5Motion(raw,p);
assert(abs(T.velocity1(end)-0.5)<1e-12);
assert(abs(T.velocity2(end)+0.3)<1e-12);
assert(abs(T.speed2(end)-0.3)<1e-12);
assert(abs(T.displacement1(end)-0.5*(time(end)-2))<1e-12);
raw(7,4)=NaN;
T=analyzeSigma5Motion(raw,p);
assert(all(isnan(T.velocity1(7:end))));
assert(all(isfinite(T.velocity2(7:end))));
disp('Independent boundary velocities, periodic unwrap and missing-data checks passed.');
end
