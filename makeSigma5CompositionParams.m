function p = makeSigma5CompositionParams(cMean,caseNumber)
% Same physical box and reference lattice for every composition.
if nargin<1, cMean=0.3; end
if nargin<2, caseNumber=1; end
validateattributes(cMean,{'numeric'},{'scalar','real','>',0,'<',1});
p=makeSigma5ShearParams(1);
p.c_init=cMean; p.setNum=100+caseNumber;
p.referenceComposition=0.3;
fractionB=3*p.referenceComposition^2-2*p.referenceComposition^3;
p.k_ini=p.k_1A+fractionB*(p.k_1B-p.k_1A);
a=2*pi/p.k_ini;
p.Nx=256; p.Ny=1024;
p.sigmaPeriods=[10 40];
p.dx=10*sqrt(10)*a/p.Nx;
p.dy=40*sqrt(10)*a/p.Ny;
p.griddim=[p.Nx p.Ny];
p.gbWidth=0.5*a;
p.disconnection.enabled=false;
p.M_c=0.1; p.M_n=1; p.beta_n=1; p.dt=0.01;
p.noiseApp_n=false; p.noiseApp_c=false;
p.shear.enabled=true; p.shear.rate=5e-4;
p.shear.every=100; p.shear.relaxSteps=2000;
p.shear.maxStrain=0.1;
p.maxStep=p.shear.relaxSteps+20000;
p.outStep=100; p.printStep=1000;
p.velocityWindow=7;
p.dispOutput=false; p.makeMovie=false;
end
