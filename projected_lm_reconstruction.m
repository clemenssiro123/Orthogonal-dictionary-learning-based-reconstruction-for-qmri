function [u, info] = projected_lm_reconstruction(data, RF, TR, r, zeta, u0, options)
%PROJECTED_LM_RECONSTRUCTION Physical L2 damping followed by box projection.
%   U=[T1,T2,density], times in ms. Frame t measures Fourier rows zeta(t):r:end.
%   Internally:
%     x=(u-[ReferenceTimes,0])./Scales, f(x)=h^2/2*||F(u)-data||^2,
%     (Hdata+lambda*G)*d=-g, x_new=clip(x+d), G=h^2*S'*S.
%   S contains the scales, so damping penalizes PHYSICAL increments S*d.
%   There is no spatial regularizer or exterior boundary condition.
%   All pixels are optimized; BodyMask is only for background diagnostics.
%
%   Approximate GN uses sparse Cholesky; exact GN uses FFT products and PCG,
%   preconditioned by the approximate matrix. LinearTolerance checks the
%   normal equations before projection. There is no descent test/backtracking.
%   lambda=max(MinDamping,InitialDamping*DampingDecay^k,weight*residual), k=0,...
%   Optional early stopping: measured residual <= DiscrepancyFactor*NoiseLevel.
%   NoiseLevel is a noise NORM, not per-component standard deviation; [] disables.
%   INFO.converged means stationarity; discrepancy_satisfied is separate.
%   Optional MonitorFcn(u) records diagnostics without changing the iteration.
%   info.time includes setup and monitoring, as in the dictionary solver.
%   Options/defaults are below; see demo_projected_lm_reconstruction.
%
%   Dong et al. (2019), sec. 4.4: https://arxiv.org/pdf/1903.11979
%   This scheduled variant does not implement Kanzow et al.'s globalization:
%   https://www.mathematik.uni-wuerzburg.de/fileadmin/10040700/paper/CNSLMP.pdf
%
%   Inputs: data is N-by-M-by-L complex Fourier data; u0 is N-by-M-by-3
%   physical [T1,T2,density]. RF/TR/zeta each describe L acquisition frames.
%   Output u has the same shape and units as u0. info contains the initial
%   state, each completed update, linear-solve accuracy and stopping reason.
%   Read sections 1-4 for setup, section 5 for updates, and section 6 for output.

    if nargin < 7, options = struct(); end
    root = fileparts(mfilename('fullpath'));
    % Make the objective, Fourier operators and Bloch simulator available.
    addpath(fullfile(root,'helper_functions'),fullfile(root,'bloch_simulations'));

    %% 1. Read settings and check the experiment
    [settings,discrepancy_threshold] = read_options(options,data,RF,TR,r,zeta,u0);
    use_exact = strcmp(settings.GaussNewton,'exact');
    % 'exact' refers to sampled Gauss-Newton J'*J products. Both modes omit
    % residual-weighted second derivatives of the nonlinear Bloch signal.
    [n,m,~] = size(u0);
    RF = RF(:).'; TR = TR(:).'; zeta = zeta(:).'; % Preserve complex RF phases.
    background = find(~settings.BodyMask); % Used only when recording diagnostics.

    %% 2. Express maps and bounds in scaled coordinates
    timer = tic; P = n*m; % Start timing and count the pixels in one map; the unknown has 3*P entries.
    scale = reshape(settings.Scales,1,1,3); % Broadcast one scale over each of the three image channels.
    % Set and reshape the coordinate origin [T1ref,T2ref,0].
    reference = [settings.ReferenceTimes 0]; offset = reshape(reference,1,1,3);
    x = (u0-offset)./scale; % Convert physical maps to scaled deviations; negative time deviations are allowed.
    % Transform the lower bounds to x coordinates and repeat them at every pixel.
    lower = repmat(reshape((settings.LowerBounds-reference)./settings.Scales,1,1,3),n,m,1);
    % Transform the upper bounds in exactly the same way.
    upper = repmat(reshape((settings.UpperBounds-reference)./settings.Scales,1,1,3),n,m,1);
    lower = lower(:); upper = upper(:); % Stack all T1 entries, then T2, then density, matching x(:).
    % Check feasibility of the initial maps before taking any step.
    assert(all(x(:)>=lower & x(:)<=upper),'Initial physical maps must satisfy bounds.');

    %% 3. Build the physical damping metric and evaluate the data model
    % Physical increment du=S*d. This diagonal metric preserves the original
    % lambda*I damping in physical coordinates even when Scales is changed.
    % Build sparse h^2*S'*S, so d'*G*d equals h^2 times the squared physical increment.
    G = settings.H^2*kron(spdiags(settings.Scales(:).^2,0,3,3),speye(P));
    % Bind the experiment to the data-only objective; no spatial penalty is added.
    model = @(z) qmri_data_objective(z,data,RF,TR,r,zeta,settings.Scales,settings.H,settings.ReferenceTimes);
    data_product = []; % Reserve the exact GN product handle used only in exact mode.
    if use_exact
        % Compute the data objective, scaled gradient, local Gram matrix and exact GN product.
        [f,~,g,Happrox,data_product] = model(x);
    else
        [f,~,g,Happrox] = model(x); % Avoid building the exact-product cache in approximate mode.
    end
    % Stop immediately if the initial objective or gradient is invalid.
    assert(isfinite(f) && all(isfinite(g)),'Initial model or derivatives are nonfinite.');

    %% 4. Allocate diagnostics for states and completed updates
    % Histories of quantities defined at each iterate, including the initial guess.
    state_fields = {'objective','data_residual_norm','projected_gradient','time'};
    % Histories describing completed updates and their linear solves.
    step_fields = {'damping','step_norm','step_norm_absolute','projection_fraction', ...
        'linear_relative_residual','linear_iterations','linear_exitflag'};
    % Preallocate space for the initial state plus all possible updates.
    for name = state_fields, info.(name{1}) = zeros(settings.MaxIterations+1,1); end
    % Preallocate per-update diagnostics to avoid growing arrays during the loop.
    for name = step_fields, info.(name{1}) = zeros(settings.MaxIterations,1); end
    % Background RMS is undefined when no background pixels are supplied.
    info.background_rms = NaN(settings.MaxIterations+1,3);
    info.monitor = zeros(settings.MaxIterations+1,0); % Filled only when a callback is supplied.
    info.stop_reason = 'max_iterations'; % Default termination reason unless another stopping condition is reached.
    info.last_linear_exitflag = NaN; % No linear solve has been attempted yet.
    % Reserve the residual of the last attempted solve, including failures.
    info.last_linear_relative_residual = NaN;
    info.last_damping = NaN; % Reserve the damping of the last attempted update.
    info.discrepancy_threshold = discrepancy_threshold; % Return the noise threshold so the demo can display it.
    % Here accepted counts completed updates, not successful objective-decrease tests.
    accepted = 0; linear_solves = 0;
    % Initialize shared variables used by the nested Hessian and preconditioner functions.
    penalty = []; Lchol = []; permutation = [];
    record_state(); % Store the current objective, stationarity measure and background values.

    %% 5. Solve damped normal equations, then project onto physical bounds
    % A linear solve first gives an unconstrained direction in scaled variables.
    % Clipping x+d enforces bounds afterwards; this is not a bound-constrained QP.
    for iteration = 1:settings.MaxIterations
        % Test whether the current residual has already reached the noise level.
        if info.data_residual_norm(accepted+1)<=discrepancy_threshold
            % Stop fitting the noisy data; this does not require a small projected gradient.
            info.stop_reason = 'discrepancy'; break;
        elseif info.projected_gradient(accepted+1)<=settings.Tolerance % Alternatively stop when the scaled box-stationarity measure is small.
            % Stop because the box-constrained stationarity tolerance is satisfied.
            info.stop_reason = 'projected_gradient'; break;
        end
        % Choose the largest of the damping floor, geometric schedule and residual-based contribution.
        damping = max([settings.MinDamping,settings.InitialDamping*settings.DampingDecay^accepted, ...
            settings.ResidualDampingWeight*info.data_residual_norm(accepted+1)]); % Use the measured residual norm for the third damping contribution.
        info.last_damping = damping; % Keep the attempted damping even if the linear solve fails.
        penalty = damping*G; % This is only the LM increment penalty; it is not part of the data objective.
        % Build the approximate damped matrix for a direct solve or exact-mode preconditioning.
        A = Happrox+penalty;
        linear_solves = linear_solves+1; % Count this linear-solve attempt.
        % The approximate damped Hessian is SPD. In exact mode this factor
        % is used ONLY for preconditioning the sampled Fourier operator.
        % Factor A(p,p)=Lchol*Lchol'; the sparse permutation reduces factor fill.
        [Lchol,chol_flag,permutation] = chol(A,'lower','vector');
        if chol_flag~=0 % Stop if numerical Cholesky factorization cannot establish positive definiteness.
            info.last_linear_exitflag = -1; % Use -1 here to identify a failed factorization.
            info.last_linear_relative_residual = Inf;
            % Keep the last valid reconstruction instead of taking an inaccurate step.
            info.stop_reason = 'linear_solve_failed'; break;
        end
        if use_exact
            % Solve the exact sampled normal equations without assembling a dense Hessian.
            [direction,flag,~,linear_iterations] = pcg(@apply_exact_hessian,-g, ...
                settings.LinearTolerance,settings.LinearMaxIterations,@apply_preconditioner); % Supply PCG's tolerance, iteration cap and approximate-matrix inverse action.
        else
            % In approximate mode the Cholesky solves directly give the desired direction.
            direction = apply_preconditioner(-g);
            % Record success and one direct solve using the same diagnostic fields as PCG.
            flag = 0; linear_iterations = 1;
        end
        % Check the actual normal equations, including FFT sampling in exact mode.
        if use_exact
            % Check the exact system A_exact*d=-g, not merely its preconditioner.
            linear_residual = apply_exact_hessian(direction)+g;
        else
            linear_residual = A*direction+g; % Check the approximate system that was solved in this mode.
        end
        % Measure relative normal-equation accuracy BEFORE clipping the direction.
        relative_residual = norm(linear_residual)/norm(g);
        info.last_linear_exitflag = flag; % Keep the solver status even when its direction is rejected.
        info.last_linear_relative_residual = relative_residual;
        % Use actual linear-system accuracy to decide whether this direction is usable.
        if ~isfinite(relative_residual) || relative_residual>settings.LinearTolerance
            % Keep the last valid reconstruction instead of taking an inaccurate step.
            info.stop_reason = 'linear_solve_failed'; break;
        end

        proposal = x(:)+direction; % Form the unconstrained LM candidate in scaled coordinates.
        % Project onto the box by clipping each entry to its allowed interval.
        candidate = min(upper,max(lower,proposal));
        step = candidate-x(:);  % Measure the actual, projected update.
        if all(step==0) % Detect a projected update that leaves every parameter unchanged.
            % Stop without calling this stall convergence to stationarity.
            info.stop_reason = 'projection_stalled'; break;
        end
        % Measure the projected step in physical L2 units, with the h^2 grid weight.
        step_squared = max(0,step'*(G*step));
        % Normalize the actual step by one plus the current deviation norm in G.
        relative_step = sqrt(step_squared)/(1+sqrt(max(0,x(:)'*(G*x(:)))));
        trial_x = reshape(candidate,n,m,3); % Convert the stacked candidate vector back to three parameter maps.
        if use_exact
            % Evaluate the candidate and prepare its next exact GN product.
            [next_f,~,next_g,next_H,next_product] = model(trial_x);
        else
            % Evaluate the candidate without building the exact-product cache.
            [next_f,~,next_g,next_H] = model(trial_x);
        end
        % Do not replace the current iterate by one with invalid objective or derivatives.
        if ~isfinite(next_f) || any(~isfinite(next_g))
            % Return the preceding valid reconstruction if candidate evaluation fails.
            info.stop_reason = 'nonfinite_model'; break;
        end
        % Scheduled projected LM: no objective-decrease acceptance condition.
        % Commit the projected update and its model values; f is allowed to increase.
        x = trial_x; f = next_f; g = next_g; Happrox = next_H;
        % Move the exact Hessian-product cache to the new linearization point.
        if use_exact, data_product = next_product; end
        accepted = accepted+1; % Count this completed projected update.
        info.damping(accepted) = damping; % Record the damping used for this completed update.
        info.step_norm(accepted) = relative_step; % Record the relative step used by the optional small-step stop.
        info.step_norm_absolute(accepted) = sqrt(step_squared); % Record the unnormalized step norm.
        % Record the fraction of parameter entries changed by clipping.
        info.projection_fraction(accepted) = nnz(candidate~=proposal)/(3*P);
        % Store the linear solve's accuracy, before projection.
        info.linear_relative_residual(accepted) = relative_residual;
        info.linear_iterations(accepted) = linear_iterations; % Store PCG iterations, or one for the direct solve.
        info.linear_exitflag(accepted) = flag; % Store the linear solver's reported status.
        record_state(); % Store the current objective, stationarity measure and background values.
        if settings.Verbose
            fprintf('%3d  f=%.6g  residual=%.6g  PG=%.2g  step=%.2g  lambda=%.2g  linres=%.2g  clipped=%.1f%%\n', ...
                accepted,f,info.data_residual_norm(accepted+1),info.projected_gradient(accepted+1), ...
                relative_step,damping,relative_residual,100*info.projection_fraction(accepted));
        end
        % Test whether the current residual has already reached the noise level.
        if info.data_residual_norm(accepted+1)<=discrepancy_threshold
            % Stop fitting the noisy data; this does not require a small projected gradient.
            info.stop_reason = 'discrepancy'; break;
        elseif info.projected_gradient(accepted+1)<=settings.Tolerance % Alternatively stop when the scaled box-stationarity measure is small.
            % Stop because the box-constrained stationarity tolerance is satisfied.
            info.stop_reason = 'projected_gradient'; break;
        elseif settings.StepTolerance>0 && relative_step<=settings.StepTolerance % Allow optional small-step stopping; it is disabled when StepTolerance is zero.
            info.stop_reason = 'small_step'; break; % A small step stops the run but does not certify stationarity.
        end
    end

    %% 6. Return physical maps, completed histories and the stopping reason
    for name = [state_fields {'background_rms','monitor'}] % Remove unused preallocated state-history entries.
        info.(name{1}) = info.(name{1})(1:accepted+1,:); % Keep the initial state and the completed updates.
    end
    % Discard unused preallocated step-history entries.
    for name = step_fields, info.(name{1}) = info.(name{1})(1:accepted,:); end
    u = offset+x.*scale; % Convert the final scaled deviations back to physical T1, T2 and density.
    info.iterations = accepted;
    info.total_linear_solves = linear_solves;
    % Only a stationarity stop sets this flag, not a discrepancy or small-step stop.
    info.converged = strcmp(info.stop_reason,'projected_gradient');
    % Report the noise-level criterion separately from stationarity.
    info.discrepancy_satisfied = info.data_residual_norm(end)<=discrepancy_threshold;
    info.elapsed_seconds = toc(timer); info.options = settings;

    %% Helper: exact sampled Gauss-Newton action for PCG
    function products = apply_exact_hessian(directions)
        % PCG requests this operation instead of a stored full data Hessian.
        % The FFT product reshapes directions into image arrays, requiring full storage.
        directions = full(directions);
        % Apply the exact sampled data GN Hessian plus physical L2 damping.
        products = data_product(directions)+penalty*directions;
    end

    %% Helper: solve with the sparse Cholesky factors
    function result = apply_preconditioner(rhs)
        % Applying the approximate inverse accelerates PCG without changing
        % the exact system that PCG is asked to solve.
        result = zeros(size(rhs)); % Allocate the answer in the original variable ordering.
        % Permute the RHS, solve L*y=b then L'*z=y, and undo the permutation.
        result(permutation,:) = Lchol'\(Lchol\rhs(permutation,:));
    end

    %% Helper: objective, stationarity and optional diagnostics
    function record_state()
        k = accepted+1; % History index 1 is the initial state, before any update.
        info.objective(k) = f; % Store the true nonlinear objective at the current iterate.
        % Recover ||Q(u)-data|| from f=h^2/2*||residual||^2 using measured samples only.
        info.data_residual_norm(k) = sqrt(2*f)/settings.H;
        % Compute x-P_box(x-g) using Euclidean clipping in scaled coordinates.
        pg = x(:)-min(upper,max(lower,x(:)-g));
        % Largest residual component; zero expresses box stationarity, not an H1 projection.
        info.projected_gradient(k) = norm(pg,inf);
        if ~isempty(background) % Compute background diagnostics only when background pixels exist.
            maps = reshape(x,P,3); % View the three maps as columns without changing their ordering.
            % Recover physical values at the background pixels without modifying x.
            values = maps(background,:).*settings.Scales+reference;
            % Store RMS background T1, T2 and density separately.
            info.background_rms(k,:) = sqrt(mean(values.^2,1));
        end
        if ~isempty(settings.MonitorFcn)
            metrics = settings.MonitorFcn(offset+x.*scale); % Truth may be used for evaluation only.
            info.monitor(k,1:numel(metrics)) = metrics(:).';
        end
        info.time(k) = toc(timer);
    end
end

%% Settings and input checks (not part of the numerical iteration)
function [settings,discrepancy_threshold] = read_options(options,data,RF,TR,r,zeta,u0)
    % assert(condition,message) stops with that message when condition is false.
    % Default coordinate scales for T1, T2 and density; physical L2 damping accounts for them.
    defaults = struct('Scales',[5000 500 100], ...
        'ReferenceTimes',[1000 100],'H',1, ...
        'GaussNewton','approximate','InitialDamping',[],'DampingDecay',.15, ... % Select the data Hessian and geometric damping schedule; [] later becomes r^2.
        'ResidualDampingWeight',1e-15,'MinDamping',1e-8, ... % Set the residual-dependent damping contribution and its positive lower limit.
        'LinearTolerance',1e-8,'LinearMaxIterations',500, ... % Set relative linear-system accuracy and the iteration limit for exact-mode PCG.
        'MaxIterations',100,'Tolerance',1e-6,'StepTolerance',0, ... % Set the outer budget and stationarity tolerance; disable small-step stopping by default.
        'NoiseLevel',[],'DiscrepancyFactor',1.1, ... % Optional noise norm delta and multiplier tau; [] disables discrepancy stopping.
        'LowerBounds',[100 10 0],'UpperBounds',[5500 600 110], ... % Physical box bounds on [T1,T2,density].
        'BodyMask',[],'Verbose',true,'MonitorFcn',[]);
    names = fieldnames(options);
    for j = 1:numel(names)
        if ~isfield(defaults,names{j}) % Reject misspelled or unsupported option names.
            error('projected_lm_reconstruction:UnknownOption','Unknown option: %s.',names{j});
        end
        defaults.(names{j}) = options.(names{j});
    end
    settings = defaults;
    % Check which data-Hessian model to use.
    settings.GaussNewton = validatestring(settings.GaussNewton,{'approximate','exact'},mfilename,'GaussNewton');
    % Check dimensions and physical bounds; MATLAB handles ordinary type errors.
    [n,m,channels] = size(u0); L = size(data,3); % Read the image dimensions, number of maps and number of frames.
    % Require three parameter maps on the same spatial grid as the data.
    assert(channels==3 && ndims(u0)==3 && size(data,1)==n && size(data,2)==m, ...
        'Expected N-by-M-by-3 maps and N-by-M-by-L data.');
    % Require a pulse, positive repetition time and sampling offset for every frame.
    assert(numel(RF)==L && numel(TR)==L && numel(zeta)==L && all(TR(:)>0), ...
        'RF, positive TR, and zeta must have one entry per frame.');
    % Check that the periodic row-sampling pattern fits the image.
    assert(r>=1 && r==fix(r) && mod(n,r)==0 && ...
        all(zeta(:)>=1 & zeta(:)<=r & zeta(:)==fix(zeta(:))), ...
        'r must divide the row count; sampling offsets must be integers in 1:r.');
    for name = {'Scales','LowerBounds','UpperBounds'} % Normalize the shape of the three per-map option vectors.
        settings.(name{1}) = reshape(settings.(name{1}),1,3); % Store the T1, T2 and density entries as a row vector.
    end
    settings.ReferenceTimes = reshape(settings.ReferenceTimes,1,2); % Store the two relaxation-time references as a row vector.
    % Ensure valid coordinate scales and physically admissible bound intervals.
    assert(all(settings.Scales>0) && all(settings.UpperBounds>settings.LowerBounds) && ...
        all(settings.LowerBounds(1:2)>0) && settings.LowerBounds(3)>=0, ...
        'Use positive scales and T1/T2 bounds, and nonnegative density bounds.');
    % Require the reference T1/T2 values to lie inside the allowed intervals.
    assert(all(settings.ReferenceTimes>=settings.LowerBounds(1:2) & settings.ReferenceTimes<=settings.UpperBounds(1:2)), ...
        'ReferenceTimes must lie within the physical bounds.');
    % A diagnostic mask, when supplied, must match the image dimensions.
    assert(isempty(settings.BodyMask) || isequal(size(settings.BodyMask),[n m]), ...
        'BodyMask must match the image size.');

    % Use the pasted method's r^2 choice when no initial damping is specified.
    if isempty(settings.InitialDamping), settings.InitialDamping = r^2; end
    % Require a positive damped system and a nonincreasing geometric schedule.
    assert(settings.H>0 && settings.InitialDamping>0 && settings.MinDamping>0 && ...
        settings.DampingDecay>0 && settings.DampingDecay<=1 && settings.ResidualDampingWeight>=0, ...
        'Use positive grid spacing/damping, decay in (0,1], and nonnegative residual weight.');
    % The relative linear-system residual tolerance must lie between zero and one.
    assert(settings.LinearTolerance>0 && settings.LinearTolerance<1, ...
        'LinearTolerance must lie in (0,1).');
    % Allow a residual slightly larger than the estimated noise norm.
    assert(settings.DiscrepancyFactor>1,'DiscrepancyFactor must be greater than 1.');
    % Comparisons with NaN are false, disabling the noise test unless it is configured.
    discrepancy_threshold = NaN;
    if ~isempty(settings.NoiseLevel) % Enable discrepancy stopping only when a noise norm is supplied.
        % Require a usable noise norm over the measured data entries.
        assert(isscalar(settings.NoiseLevel) && isfinite(settings.NoiseLevel) && settings.NoiseLevel>=0, ...
            'NoiseLevel must be a finite, nonnegative noise norm.');
        % Compute tau*delta in the same unweighted norm as data_residual_norm.
        discrepancy_threshold = settings.DiscrepancyFactor*settings.NoiseLevel;
        % Reject overflow in the stopping threshold.
        assert(isfinite(discrepancy_threshold),'Discrepancy threshold must be finite.');
    end

end
