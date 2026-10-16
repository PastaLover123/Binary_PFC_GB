function simulateXPFC(params)
% Run with defaults, or supply a complete parameter structure.
if nargin < 1
    params = genParams();
end
if ~isfolder(params.output_dir), mkdir(params.output_dir); end
%% Initialize GPU
tic;
disp('[[[[[ Preparing for simulation run... ]]]]]');
disp(['[[[[[ Run ID is ' num2str(params.setNum) ' ]]]]]']);
toc;
%% Prepare operators on GPU
tic;
disp('[[[[[ Preparing operators... ]]]]]');
[ operators] = genOperators( params );
toc;
%% Prepare output files
tic;
disp('[[[[[ Preparing output files... ]]]]]');
filename_params = strcat('o_',num2str(params.setNum),'_params.mat');
fullfilename_params = fullfile(params.output_dir,filename_params);
save(fullfilename_params,'params');

if params.makeMovie
    filename_video1=strcat('s',num2str(params.setNum),'_simulmov_n.avi');
    fullfilename_video1 = fullfile(params.output_dir,filename_video1);
    v1 = VideoWriter(fullfilename_video1);
    v1.FrameRate = params.movieFramerate;
    open(v1)
    filename_video2=strcat('s',num2str(params.setNum),'_simulmov_c.avi');
    fullfilename_video2 = fullfile(params.output_dir,filename_video2);
    v2 = VideoWriter(fullfilename_video2);
    v2.FrameRate = params.movieFramerate;
    open(v2)
end
toc;
%% Prepare simulation fields
tic;
disp('[[[[[ Preparing simulation fields... ]]]]]');
fields = struct();
% Save n, c fields and their Fourier transforms to gpuArrays
fields.n = zeros(params.griddim,"gpuArray");
fields.c = zeros(params.griddim,"gpuArray");
fields.n_k = complex(zeros(params.griddim,"gpuArray"));
fields.c_k = complex(zeros(params.griddim,"gpuArray"));
% Keep track of simulation parameters. If any of these are of order system
% size rather than scalars, they should be made gpuArrays.
status = struct();
status.timestep = 0;
status.debug = 0;
toc;
%% Initialize simulation fields
tic;
disp('[[[[[ Initializing simulation fields... ]]]]]');

% Prepare the density and composition fields in real-space
[fields] = genInitialConditions(fields,params);

fields.n_k = fftn(fields.n);
fields.c_k = fftn(fields.c);
if isfield(operators,'spectralMask')
    discarded = fields.n_k.*(1-operators.spectralMask);
    fraction = gather(sum(abs(discarded(:)).^2) ...
        /max(sum(abs(fields.n_k(:)).^2),eps));
    fprintf('Initial density spectral power removed at Nyquist: %.3e\n',fraction);
    fields.n_k = fields.n_k.*operators.spectralMask;
    fields.c_k = fields.c_k.*operators.spectralMask;
    fields.n = real(ifftn(fields.n_k));
    fields.c = real(ifftn(fields.c_k));
end
% Equal density histories initialize zero discrete density velocity.
fields.n_k_prev = fields.n_k;

if params.dispOutput
    figure(1)
    subplot(1,2,1)
    imagesc(fields.n);
    colorbar;
    pause(.1);
    subplot(1,2,2)
    imagesc(fields.c);
    colorbar;
    pause(.1);
end
toc;
trackGB = isfield(params,'trackGB') && params.trackGB;
gbRows = zeros(0,8);
lastWrapped = []; lastUnwrapped = [];
if trackGB
    [pos,profile,indicator,localGB] = measureSigma5(fields,params);
    status.localGB = localGB;
    status.localGBx = ((0:params.Nx-1)'-floor(params.Nx/2))*params.dx;
    if params.dispOutput
        figure(4);
        for boundary=1:2
            subplot(2,1,boundary);
            plot(status.localGBx,localGB(:,boundary));
            xlabel('Material X'); ylabel(sprintf('GB %d: Y(X)',boundary));
        end
        drawnow;
    end
    if any(~isfinite(pos)), error('Initial Sigma5 boundaries not resolved.'); end
    initialPositions = pos;
    lastWrapped = pos; lastUnwrapped = pos;
    gbRows = [0 params.gamma pos 0 0 NaN NaN];
end
fieldRows = diagnosticRow(fields,params,0);
if isfield(params,'saveFields') && params.saveFields
    saveState(fields,params,status);
end
%% Evolve fields
tic;
disp('[[[[[ Starting main simulation loop... ]]]]]');
for step=1:params.maxStep
    status.timestep = step;
    if isfield(params,'freezeCompositionSteps') ...
            && step == params.freezeCompositionSteps+1
        fprintf('Composition released at step %d (time %.6g).\n',step,(step-1)*params.dt);
    end
    if isfield(params,'shear') && params.shear.enabled ...
            && step == params.shear.relaxSteps+1
        fprintf('Starting shear stage at step %d.\n',step);
    end
    % Piecewise affine shear of the periodic cell; fields stay on material grid.
    if isfield(params,'shear') && params.shear.enabled ...
            && step > params.shear.relaxSteps ...
            && mod(step-params.shear.relaxSteps,params.shear.every)==0
        increment = params.shear.rate*params.dt*params.shear.every;
        gammaNew = max(-params.shear.maxStrain, ...
            min(params.shear.maxStrain, params.gamma+increment));
        if gammaNew ~= params.gamma
            params.gamma = gammaNew;
            operators = genOperators(params);
        end
    end
    % Reject invalid input before evaluating composition logarithms.
    if any(~isfinite(fields.n(:))) || any(~isfinite(fields.c(:))) ...
            || any(fields.c(:) <= 0 | fields.c(:) >= 1)
        error('Invalid fields before step %d; reduce dt and restart.', step);
    end
    % Perform one time-step
    [fields,status] = stepXPFC( params, fields, operators, status);
    % Save to videos,
    if mod(step,params.outStep)==0
        % Basic error check: Look for divergences in n or c and check c is
        % within bounds ( 0 <= c <= 1)
        if any(~isfinite(fields.n(:))) || any(~isfinite(fields.c(:)))
            error(['ERROR: OBTAINED NAN VALUES BY STEP ' num2str(step)])
        end
        if nnz((fields.c <= 0) )>0 || nnz((fields.c >= 1)) > 0
            error(['ERROR: OBTAINED C FIELD OUT OF BOUNDS BY STEP ' num2str(step)])
        end
        if trackGB
            [pos,profile,indicator,localGB] = measureSigma5(fields,params);
            status.localGB = localGB;
            status.localGBx = ((0:params.Nx-1)'-floor(params.Nx/2))*params.dx;
            displacement = [NaN NaN]; velocity = [NaN NaN];
            if all(isfinite(pos)) && ~isempty(lastWrapped)
                Ly = params.Ny*params.dy;
                unwrapped = lastUnwrapped + mod(pos-lastWrapped+Ly/2,Ly)-Ly/2;
                displacement = unwrapped-initialPositions;
                lastWrapped = pos; lastUnwrapped = unwrapped;
            elseif all(isfinite(pos))
                warning('Boundary identities lost; displacements remain undefined.');
            else
                warning('Expected two GB crossings at step %d; tracking invalid.',step);
                lastWrapped = []; lastUnwrapped = [];
            end
            row = [step*params.dt params.gamma pos displacement velocity];
            gbRows(end+1,:) = row;
            % Separate signed velocities from trailing fits; never average GBs.
            if size(gbRows,1) >= 3
                ii = max(1,size(gbRows,1)-4):size(gbRows,1);
                if all(isfinite(gbRows(ii,5:6)), 'all')
                    for b=1:2
                        fit = polyfit(gbRows(ii,1),gbRows(ii,4+b),1);
                        gbRows(end,6+b) = fit(1);
                    end
                end
            end
            T = array2table([gbRows abs(gbRows(:,7:8))], 'VariableNames', ...
                {'time','gamma','y1','y2','displacement1','displacement2', ...
                 'velocity1','velocity2','speed1','speed2'});
            writetable(T,fullfile(params.output_dir,'gb_tracking.csv'));
            save(fullfile(params.output_dir,'gb_tracking.mat'),'T','profile','localGB','params');
            fprintf('GB displacements: [%g, %g]; gamma=%g\n', displacement,params.gamma);
            if params.dispOutput
                figure(3);
                subplot(2,1,1);
                plot(T.time,T.displacement1,T.time,T.displacement2);
                xlabel('Time'); ylabel('GB displacement');
                legend('GB 1','GB 2');
                subplot(2,1,2);
                plot(T.time,T.speed1,T.time,T.speed2);
                xlabel('Time'); ylabel('|GB velocity|');
                legend('GB 1','GB 2'); drawnow;
                figure(4);
                for boundary = 1:2
                    subplot(2,1,boundary);
                    plot(status.localGBx,localGB(:,boundary));
                    xlabel('Material X'); ylabel(sprintf('GB %d: Y(X)',boundary));
                end
                drawnow;
            end
        end
        fieldRows(end+1,:) = diagnosticRow(fields,params,step);
        fieldTable = array2table(fieldRows,'VariableNames', ...
            {'time','gamma','meanN','meanC','minN','maxN','minC','maxC','stdC'});
        writetable(fieldTable,fullfile(params.output_dir,'field_diagnostics.csv'));
        fprintf('Composition: mean=%.8g, range=[%.6g, %.6g], std=%.4g\n', ...
            fieldRows(end,4),fieldRows(end,7:9));
        if isfield(params,'saveFields') && params.saveFields
            saveState(fields,params,status);
        end
        % Save n and c fields to movies
        if params.makeMovie
            data1 = fields.n;
            data2 = fields.c;
            % Make the data grayscale
            minv1 = min(data1(:));
            maxv1 = max(data1(:));
            frame=gather((data1-minv1)/max(maxv1-minv1,eps));
            writeVideo(v1,frame);
            minv2 = min(data2(:));
            maxv2 = max(data2(:));
            frame=gather((data2-minv2)/max(maxv2-minv2,eps));
            writeVideo(v2,frame);
        end
        % Plot n and c fields
        if params.dispOutput
            figure(2)
            subplot(1,3,1)
            imagesc(fields.n);
            colorbar;
            pause(.1);
            subplot(1,3,2)
            imagesc(fields.c);
            colorbar;
            pause(.1);
            % Plot the strain with reference wave-vector equal to k_1A
            subplot(1,3,3)
            if trackGB
                imagesc(indicator); title('Orientation: grain 1 minus grain 2');
            else
                imagesc(real(strainFilter(sqrt(operators.k2),operators.k2,params.k_1A,1,1,1,fields.n)))
            end
            colorbar;
            pause(.1);
        end
    end
    % Output message
    if mod(step,params.printStep)==0
        timed = toc;
        etatime = (params.maxStep - step)/params.printStep*timed;
        disp(  '######################################################################');
        disp(  ['### Step ' num2str(step) ' / ' num2str(params.maxStep) ' (' num2str(round(step/params.maxStep*100)) '%) ###']  );
        disp( ['### Time since last print out: ' num2str(timed) ' seconds. ### ETA: ' char(duration(0,0,etatime)) ' ###' ] )
        tic
    end
end

if params.makeMovie
    close(v1)
    close(v2)
end

disp('[[[[[ Done. ]]]]]');
end


function row = diagnosticRow(fields,params,step)
gamma = 0;
if isfield(params,'gamma'), gamma = params.gamma; end
row = [step*params.dt gamma gather(mean(fields.n(:))) ...
    gather(mean(fields.c(:))) gather(min(fields.n(:))) ...
    gather(max(fields.n(:))) gather(min(fields.c(:))) ...
    gather(max(fields.c(:))) gather(std(fields.c(:)))];
end

function saveState(fields,params,status)
% Preserve doubles and the density history without image downsampling.
snapshot = struct();
names = fieldnames(fields);
for j = 1:numel(names)
    snapshot.(names{j}) = gather(fields.(names{j}));
end
save(fullfile(params.output_dir,sprintf('s%d_step%07d.mat', ...
    params.setNum,status.timestep)), 'snapshot','params','status','-v7.3');
end
