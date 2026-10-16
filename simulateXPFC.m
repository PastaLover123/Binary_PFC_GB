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
%% Evolve fields
tic;
disp('[[[[[ Starting main simulation loop... ]]]]]');
for step=1:params.maxStep
    status.timestep = step;
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
        if nnz(isnan(fields.n))>0||nnz(isnan(fields.c))>0
            error(['ERROR: OBTAINED NAN VALUES BY STEP ' num2str(step)])
        end
        if nnz((fields.c <= 0) )>0 || nnz((fields.c >= 1)) > 0
            error(['ERROR: OBTAINED C FIELD OUT OF BOUNDS BY STEP ' num2str(step)])
        end
        % Lossless full-grid snapshot, including density history for restarts.
        if isfield(params, 'saveFields') && params.saveFields
            snapshot = struct();
            names = fieldnames(fields);
            for j = 1:numel(names)
                snapshot.(names{j}) = gather(fields.(names{j}));
            end
            save(fullfile(params.output_dir, sprintf('s%d_step%07d.mat', ...
                params.setNum, step)), 'snapshot', 'params', 'status', '-v7.3');
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
            imagesc(real(strainFilter(sqrt(operators.k2),operators.k2,params.k_1A,1,1,1,fields.n)))
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
