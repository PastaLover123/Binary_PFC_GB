function outputDir = runSigma5CompositionSweep(cList,options)
% Sequential GPU sweep of plain Sigma5 bicrystals. No disconnections.
% outputDir=runSigma5CompositionSweep;
if nargin<1 || isempty(cList), cList=[0.1 0.2 0.3 0.4 0.5]; end
if nargin<2, options=struct(); end
validateattributes(cList,{'numeric'},{'vector','real','finite','>',0,'<',1});
if ~isfield(options,'relaxSteps'), options.relaxSteps=2000; end
if ~isfield(options,'loadingSteps'), options.loadingSteps=20000; end
if ~isfield(options,'saveSnapshots'), options.saveSnapshots=true; end
if ~isfield(options,'outputDir')
    [~,tag]=fileparts(tempname);
    options.outputDir=fullfile(pwd,['sigma5_composition_sweep_' tag]);
end
validateattributes(options.relaxSteps,{'numeric'},{'scalar','integer','nonnegative'});
validateattributes(options.loadingSteps,{'numeric'},{'scalar','integer','positive'});
assert(mod(options.relaxSteps,100)==0 && mod(options.loadingSteps,100)==0, ...
    'Relaxation and loading durations must be multiples of 100 steps.');
outputDir=options.outputDir;
if ~isfolder(outputDir), mkdir(outputDir); end
cases=cell(numel(cList),1);
fprintf('Five-composition default: 256 along X (GB tangent), 1024 along Y.\n');
for caseNumber=1:numel(cList)
    p=makeSigma5CompositionParams(cList(caseNumber),caseNumber);
    p.shear.relaxSteps=options.relaxSteps;
    p.maxStep=options.relaxSteps+options.loadingSteps;
    p.output_dir=fullfile(outputDir,sprintf('case%02d_c%.3f',caseNumber,p.c_init));
    if ~isfolder(p.output_dir), mkdir(p.output_dir); end
    caseResult=struct('composition',p.c_init,'params',p,'status','running', ...
        'message','','tracking',table());
    fprintf('\nCase %d/%d: c=%.3f; output %s\n', ...
        caseNumber,numel(cList),p.c_init,p.output_dir);
    raw=zeros(0,8); profiles=zeros(p.Ny,0);
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
