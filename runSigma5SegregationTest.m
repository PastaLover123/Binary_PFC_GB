function outputDir = runSigma5SegregationTest(cMean,options)
% One GPU Sigma5 test with stronger mismatch and narrower correlation peaks.
% folder=runSigma5SegregationTest(0.3);
% Baseline: runSigma5SegregationTest(0.3,struct("mismatch",0.05,"widthScale",1));
if nargin<1 || isempty(cMean), cMean=0.3; end
cList=cMean;
if nargin<2, options=struct(); end
validateattributes(cList,{'numeric'},{'scalar','real','finite','>',0,'<',1});
if ~isfield(options,'mismatch'), options.mismatch=0.10; end
if ~isfield(options,'widthScale'), options.widthScale=0.8; end
validateattributes(options.mismatch,{'numeric'},{'scalar','real','finite','>=',0,'<',0.3});
validateattributes(options.widthScale,{'numeric'},{'scalar','real','finite','positive'});
if ~isfield(options,'relaxSteps'), options.relaxSteps=2000; end
if ~isfield(options,'loadingSteps'), options.loadingSteps=20000; end
if ~isfield(options,'saveSnapshots'), options.saveSnapshots=true; end
if ~isfield(options,'outputDir')
    [~,tag]=fileparts(tempname);
    options.outputDir=fullfile(pwd,['sigma5_segregation_test_' tag]);
end
validateattributes(options.relaxSteps,{'numeric'},{'scalar','integer','nonnegative'});
validateattributes(options.loadingSteps,{'numeric'},{'scalar','integer','positive'});
assert(mod(options.relaxSteps,100)==0 && mod(options.loadingSteps,100)==0, ...
    'Relaxation and loading durations must be multiples of 100 steps.');
outputDir=options.outputDir;
if ~isfolder(outputDir), mkdir(outputDir); end
cases=cell(numel(cList),1);
fprintf('Single composition: 256 along X, 1024 along Y; mismatch %.3f, width scale %.3f.\n',options.mismatch,options.widthScale);
for caseNumber=1:numel(cList)
    p=makeTestParams(cList(caseNumber),options);
    p.shear.relaxSteps=options.relaxSteps;
    p.maxStep=options.relaxSteps+options.loadingSteps;
    p.output_dir=fullfile(outputDir,sprintf('case%02d_c%.3f',caseNumber,p.c_init));
    if ~isfolder(p.output_dir), mkdir(p.output_dir); end
    caseResult=struct('composition',p.c_init,'params',p,'status','running', ...
        'message','','tracking',table());
    fprintf('\nCase %d/%d: c=%.3f; output %s\n', ...
        caseNumber,numel(cList),p.c_init,p.output_dir);
    raw=zeros(0,8); profiles=zeros(p.Ny,0); segregation=zeros(0,8);
    compositionProfiles=zeros(p.Ny,0);
    try
        operators=genOperators(p);
        fields=struct('n',zeros(p.griddim,'gpuArray'),'c',zeros(p.griddim,'gpuArray'));
        fields=genInitialConditions(fields,p);
        fields.n_k=fftn(fields.n).*operators.spectralMask;
        fields.c_k=fftn(fields.c).*operators.spectralMask;
        fields.n=real(ifftn(fields.n_k)); fields.c=real(ifftn(fields.c_k));
        fields.n_k_prev=fields.n_k;
        status=struct('timestep',0);
        m=measureSigma5Profile(fields,p);
        raw(end+1,:)=sampleRow(fields,p,0,m);
        profiles(:,end+1)=m.profile;
        [segregation,compositionProfiles]=recordSegregation(fields,p,m,raw(end,1),segregation,compositionProfiles);
        if options.saveSnapshots, saveSnapshot(fields,p,0); end
        if any(~isfinite(m.positions)), error('Initial GB profile is unresolved.'); end
        timer=tic;
        for step=1:p.maxStep
            status.timestep=step;
            if any(~isfinite(fields.n(:))) || any(~isfinite(fields.c(:))) ...
                    || any(fields.c(:)<=0 | fields.c(:)>=1)
                error('Invalid fields before step %d; reduce dt and repeat this case.',step);
            end
            if step>p.shear.relaxSteps && mod(step-p.shear.relaxSteps,p.shear.every)==0
                newGamma=min(p.shear.maxStrain,p.gamma+p.shear.rate*p.dt*p.shear.every);
                if newGamma~=p.gamma
                    p.gamma=newGamma;
                    operators=genOperators(p);
                end
            end
            [fields,status]=stepXPFC(p,fields,operators,status);
            if mod(step,p.outStep)==0
                m=measureSigma5Profile(fields,p);
                raw(end+1,:)=sampleRow(fields,p,step,m);
                profiles(:,end+1)=m.profile;
                [segregation,compositionProfiles]=recordSegregation(fields,p,m,step,segregation,compositionProfiles);
                if any(~isfinite(raw(end,6:7))) ...
                        || any(~isfinite(fields.n(:))) || any(~isfinite(fields.c(:))) ...
                        || any(fields.c(:)<=0 | fields.c(:)>=1)
                    error('Invalid fields at step %d.',step);
                end
                T=analyzeSigma5Motion(raw,p);
                save(fullfile(p.output_dir,'tracking.mat'),'T','raw','profiles','p');
                writetable(T,fullfile(p.output_dir,'gb_tracking.csv'));
                caseResult.tracking=T; caseResult.params=p;
                cases{caseNumber}=caseResult;
                save(fullfile(outputDir,'composition_sweep.mat'),'cases','cList','options');
                if options.saveSnapshots && (step==p.shear.relaxSteps || step==p.maxStep)
                    saveSnapshot(fields,p,step);
                end
                if all(isfinite(m.positions))
                    width=mod(m.positions(2)-m.positions(1),p.Ny*p.dy);
                    if min(width,p.Ny*p.dy-width)<6*(2*pi/p.k_ini)
                        caseResult.status='stopped: GBs close';
                        caseResult.message=sprintf('Grain width below six lattice units at step %d.',step);
                        if options.saveSnapshots, saveSnapshot(fields,p,step); end
                        break;
                    end
                else
                    caseResult.status='stopped: tracking unresolved';
                    caseResult.message=sprintf('No unique pair of GB crossings at step %d.',step);
                    if options.saveSnapshots, saveSnapshot(fields,p,step); end
                    break;
                end
            end
            if mod(step,p.printStep)==0
                fprintf('c=%.2f: %d/%d, gamma=%.4f, elapsed %.1fs\n', ...
                    p.c_init,step,p.maxStep,p.gamma,toc(timer));
            end
        end
        if strcmp(caseResult.status,'running'), caseResult.status='completed'; end
    catch err
        caseResult.status='failed';
        caseResult.message=err.message;
        warning('Composition %.3f: %s',p.c_init,err.message);
    end
    caseResult.params=p;
    if ~isempty(raw), caseResult.tracking=analyzeSigma5Motion(raw,p); end
    cases{caseNumber}=caseResult;
    save(fullfile(outputDir,'composition_sweep.mat'),'cases','cList','options');
    fprintf('Case c=%.3f: %s\n',p.c_init,caseResult.status);
end
plotSigma5CompositionSweep(outputDir);
fprintf('Results and comparison figure: %s\n',outputDir);
end

function row=sampleRow(fields,p,step,m)
row=[step step*p.dt p.gamma m.positions ...
    gather(mean(fields.n(:))) gather(mean(fields.c(:))) m.strength];
end

function saveSnapshot(fields,p,step)
snapshot=struct();
names=fieldnames(fields);
for j=1:numel(names), snapshot.(names{j})=gather(fields.(names{j})); end
save(fullfile(p.output_dir,sprintf('fields_step%07d.mat',step)), ...
    'snapshot','p','-v7.3');
end

function p=makeTestParams(cMean,options)
p=makeSigma5CompositionParams(cMean,1);
p.a2=p.a1*(1-options.mismatch);
p.k_1A=2*pi/p.a1; p.k_1B=2*pi/p.a2;
p.k_2A=sqrt(2)*p.k_1A; p.k_2B=sqrt(2)*p.k_1B;
p.sigma_1A=p.sigma_1A*options.widthScale;
p.sigma_1B=p.sigma_1B*options.widthScale;
p.sigma_2A=sqrt(2)*p.sigma_1A;
p.sigma_2B=sqrt(2)*p.sigma_1B;
% Use the same reference-composition rule as the sweep, updating the box too.
fB=3*p.referenceComposition^2-2*p.referenceComposition^3;
p.k_ini=p.k_1A+fB*(p.k_1B-p.k_1A);
a=2*pi/p.k_ini;
p.dx=p.sigmaPeriods(1)*sqrt(10)*a/p.Nx;
p.dy=p.sigmaPeriods(2)*sqrt(10)*a/p.Ny;
p.gbWidth=0.5*a;
end

function [history,profiles]=recordSegregation(fields,p,m,step,history,profiles)
c=gather(fields.c); n=gather(fields.n); profile=mean(c,1).';
profiles(:,end+1)=profile;
y=((0:p.Ny-1)'-floor(p.Ny/2))*p.dy;
Ly=p.Ny*p.dy; a=2*pi/p.k_ini;
distance=mod(y-m.positions+Ly/2,Ly)-Ly/2;
bulk=all(abs(distance)>5*a,2);
cBulk=NaN; excess=[NaN NaN]; peak=[NaN NaN];
if all(isfinite(m.positions)) && any(bulk)
    cBulk=mean(profile(bulk));
    for b=1:2
        window=abs(distance(:,b))<=3*a;
        excess(b)=sum(profile(window)-cBulk)*p.dy;
        peak(b)=max(profile(window));
    end
end
history(end+1,:)=[step step*p.dt cBulk excess peak mean(c(:))];
S=array2table(history,'VariableNames',{'step','time','bulkC','excess1','excess2','peakMeanProfile1','peakMeanProfile2','meanC'});
writetable(S,fullfile(p.output_dir,'segregation.csv'));
save(fullfile(p.output_dir,'segregation.mat'),'S','profiles','y','p');
fig=findall(groot,'Type','figure','Tag','sigma5Segregation');
if isempty(fig), fig=figure('Tag','sigma5Segregation','Name','Sigma5 segregation'); else, fig=fig(1); end
% Keep the initial density scale: autoscaling would hide melting.
if step==0 || ~isappdata(fig,'densityLimits')
    setappdata(fig,'densityLimits',[min(n(:))-eps max(n(:))+eps]);
    setappdata(fig,'lastPositions',[-Ly/4 Ly/4]);
end
if all(isfinite(m.positions)), setappdata(fig,'lastPositions',m.positions); end
plotPositions=getappdata(fig,'lastPositions');
plotDistance=mod(y-plotPositions+Ly/2,Ly)-Ly/2;
figure(fig); clf(fig); set(fig,'Position',[100 100 1400 1000]);
x=(0:p.Nx-1)*p.dx;
for b=1:2
    [offset,order]=sort(plotDistance(:,b)); keep=abs(offset)<=4*a;
    for channel=1:2
        ax=subplot(3,2,2*(b-1)+channel);
        if channel==1
            values=n; limits=getappdata(fig,'densityLimits'); label='density n';
            colormap(ax,gray(256));
        else
            values=c; limits=[min(c(:))-eps max(c(:))+eps]; label='composition c';
            colormap(ax,parula(256));
        end
        imagesc(x,offset(keep)/a,values(:,order(keep)).'); axis xy;
        caxis(limits); colorbar;
        xlabel('X'); ylabel('(Y - Y_{GB}) / a');
        title(sprintf('GB %d: %s, step %d',b,label,step));
        if any(~isfinite(m.positions))
            title(sprintf('GB %d: %s (last tracked location), step %d',b,label,step));
        end
    end
end
subplot(3,2,5); hold on;
for b=1:2
    [offset,order]=sort(distance(:,b)); keep=abs(offset)<=4*a;
    plot(offset(keep)/a,profile(order(keep)),'DisplayName',sprintf('GB %d',b));
end
yline(p.c_init,':','Initial mean'); xlabel('(Y - Y_{GB}) / a'); ylabel('Mean composition along X'); legend('Location','best');
subplot(3,2,6); plot(history(:,2),history(:,4:5));
xlabel('Time'); ylabel('Integrated excess composition'); legend('GB 1','GB 2','Location','best');
drawnow;
savefig(fig,fullfile(p.output_dir,'segregation.fig'));
exportgraphics(fig,fullfile(p.output_dir,'segregation.png'),'Resolution',300);
end
