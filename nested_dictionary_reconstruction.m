function [u,D,C,info] = nested_dictionary_reconstruction(data,RF,TR,r,zeta,u0,options)
%NESTED_DICTIONARY_RECONSTRUCTION Alternate dictionary learning and map updates.
%   u0 and u contain [T1,T2,density]; relaxation times are in milliseconds.
%   data contains Fourier measurements; RF/TR describe the pulse sequence.
%   Frame t measures rows zeta(t):r:end. options holds the solver settings.
%
%   Read the numbered sections below in this order:
%     1-5: settings, scaled maps, patch geometry, initial D/C and QP setup.
%       6: repeat A) dictionary learning, B) bounded map update,
%                 C) accept the complete update, D) check stopping.
%       7: convert to physical maps and return the recorded histories.
%   Local helpers below the algorithm evaluate models, record diagnostics,
%   and check settings. No ground-truth maps are used by the algorithm.
%
%   The objective, with x=(u-[ReferenceTimes,0])./Scales, is
%
%     h^2/2*||Q(u)-data||^2 + Alpha/2*x'*K*x
%       + sum_j Lambda(j)*(0.5*||R*x_j-D_j*C_j||_F^2 + Beta(j)*||C_j||_1).
%
%   R extracts overlapping patches; each column of C represents one patch.
%   Each D(:,:,j) is orthogonal. K=h^2*diag(-LAP,-LAP,-LAP) penalizes spatial
%   gradients, with zero exterior values for x. Beta is INSIDE Lambda*(...).
%   This is the discrete P2/Appendix A convention used by this project.
%
%   InnerStopRule: 'tolerance' = nested; 'one-step' = one D/C pair per map;
%   'fixed' = a prescribed number of pairs. The map update is the same.
%   The outer stop requires BOTH small map and small total D/C changes.
%   A failed inner solve or map update returns the last COMPLETE iterate.
%   BodyMask/MonitorFcn are diagnostics only; all pixels remain variables.
%   MATLAB R2021b + Optimization Toolbox. See demo_compare_all_reconstructions.
%
%   Array sizes: data is N-by-M-by-L; u/u0 are N-by-M-by-3 physical maps.
%   For A=prod(PatchSize) and Q overlapping patches, D is A-by-A-by-3 and
%   C is A-by-Q-by-3. Channel order is always T1,T2,rho. Returned D/C describe
%   patches of scaled deviations x, not physical u. State histories include
%   the initial point; step histories describe only completed outer updates.

    if nargin<7, options=struct(); end
    root=fileparts(mfilename('fullpath'));
    addpath(fullfile(root,'helper_functions'),fullfile(root,'bloch_simulations'));

    %% 1. Settings: choose nested or one-step dictionary learning
    settings=read_options(options,data,RF,TR,r,zeta,u0);
    inner_rule=settings.InnerStopRule;
    inner_limit=settings.InnerMaxIterations;
    if strcmp(settings.InnerStopRule,'one-step')
        inner_rule='fixed'; % The inner solver uses fixed-count mode for one-step DL.
        inner_limit=1;      % Exactly one D update and one C update per map.
    end
    use_exact=strcmp(settings.GaussNewton,'exact');
    % This changes data curvature only. Dictionary updates and the nonlinear
    % objective/gradient use the same formulas in both GN modes.
    [n,m,~]=size(u0);
    P=n*m; % Pixels per parameter map; the map QP has at most 3*P unknowns.
    RF=RF(:).'; TR=TR(:).'; zeta=zeta(:).'; % .' preserves complex RF phases.

    %% 2. Scale physical maps and translate the physical bounds
    timer=tic;
    scale=reshape(settings.Scales,1,1,3); % One scale per channel, broadcast over pixels.
    offset=reshape([settings.ReferenceTimes 0],1,1,3); % Origin [T1ref,T2ref,0].
    x=(u0-offset)./scale; % Dimensionless deviations; these may be negative.
    % Bounds must be expressed in the same coordinates as the QP unknowns.
    lower=repmat((reshape(settings.LowerBounds,1,1,3)-offset)./scale,n,m,1);
    upper=repmat((reshape(settings.UpperBounds,1,1,3)-offset)./scale,n,m,1);
    lower=lower(:); upper=upper(:); % Map-major ordering, identical to x(:).
    % Stop here if the starting maps are nonfinite or physically infeasible.
    assert(all(isfinite(x(:)) & x(:)>=lower & x(:)<=upper),'Initial maps must satisfy bounds.');
    free=lower<upper; % Equal-bound variables have zero increments, not zero physical values.
    number_of_free=nnz(free); % Fixed-bound variables never enter quadprog.

    %% 3. Build patch extraction, spatial regularization and the damping metric
    % R*x(:) stacks valid patches. counts tells how many patches cover each pixel.
    [R,number_of_patches,counts]=create_R_matrix([n m],settings.PatchSize,1);
    patch_area=prod(settings.PatchSize); % Number of entries in a single patch.
    % x'*K*x is the summed gradient energy of all three scaled maps.
    K=settings.H^2*kron(speye(3),-dirichlet_laplacian_matrix(n,m,settings.H));
    G=K; % H01 damping measures spatial gradients of the map increment.
    % Full H1 damping also measures the increment values themselves.
    if strcmp(settings.StepNorm,'H1'), G=G+settings.H^2*speye(3*P); end
    patch_diagonal=counts*settings.Lambda; % One overlap-count vector per parameter.
    % Outer product: entry (p,j) is counts(p)*Lambda(j), a P-by-3 array.
    % R'*R is diagonal, so the patch-fit Hessian needs only overlap counts.
    regularizer_H=settings.Alpha*K+spdiags(patch_diagonal(:),0,3*P,3*P);

    %% 4. Initialize an orthogonal dictionary and sparse coefficients per map
    D=settings.D0;
    C=settings.C0;
    if isempty(D), D=zeros(patch_area,patch_area,3); end
    if isempty(C), C=zeros(patch_area,number_of_patches,3); end
    assert(isequal(size(D),[patch_area patch_area 3]) && ...
        isequal(size(C),[patch_area number_of_patches 3]), 'D0 or C0 has incompatible size.');
    for j=1:3
        X=create_data_matrix_fast(x(:,:,j),settings.PatchSize,R); % One patch per column.
        % PCA of the starting patches supplies an orthogonal basis if D0 is absent.
        if isempty(settings.D0), [D(:,:,j),~,~]=svd(X*X','econ'); end
        % Orthogonality is required for the closed-form D/C update formulas.
        assert(norm(D(:,:,j)'*D(:,:,j)-eye(patch_area),'fro')<1e-8*sqrt(patch_area), ...
            'Each initial dictionary must be orthogonal.');
        % Project onto the basis and soft-threshold to initialize sparse codes.
        if isempty(settings.C0), C(:,:,j)=prox_c(settings.Beta(j),D(:,:,j)'*X,'soft'); end
    end
    initial_dictionary_norm=hypot(norm(D(:)),norm(C(:))); % Fixed scale in eta_k.
    % Keep this scale fixed even when the norms of later D/C iterates change.

    %% 5. Prepare data derivatives, the bounded QP solver and output histories
    % This function evaluates the sampled Fourier/Bloch data-fit term at maps z.
    data_model=@(z) qmri_data_objective(z,data,RF,TR,r,zeta,settings.Scales,settings.H,settings.ReferenceTimes);
    data_product=[]; penalty=[]; % Shared by the exact QP Hessian callback.
    fdata=[]; gdata=[]; Happrox=[]; % Values updated by the nested data-model function.
    targets=zeros(size(C)); % Cache the reconstructed patches D*C.
    sparsity_value=0;       % Cache sum_j Lambda(j)*Beta(j)*||C_j||_1.
    update_patch_targets(D,C); % Cache D*C and the coefficient penalty for fixed dictionaries.
    update_data_model(); % Data derivatives depend on x, not on D or C.
    qpopts=optimoptions('quadprog','Algorithm','trust-region-reflective','Display','off', ...
        'OptimalityTolerance',settings.QPOptimalityTolerance,'FunctionTolerance',settings.QPFunctionTolerance, ...
        'TolPCG',settings.QPPCGTolerance,'PrecondBandWidth',inf,'MaxIterations',settings.QPMaxIterations); 
    % The supplied sparse Hessian also preconditions exact FFT products.

    % Exact mode uses FFT products for the data Hessian; the sparse matrix preconditions it.
    if use_exact, qpopts.HessianMultiplyFcn=@exact_hessian; end

    [info,state_fields,step_fields]=initialize_history(settings.MaxIterations);
    completed=0;  % Count only complete, accepted outer updates.
    total_trials=0; % Count every map QP attempt, including rejected attempts.
    damping=max(settings.InitialDamping,settings.MinDamping);
    record_state(); % History entry 1 describes the supplied starting maps.

    %% 6. Alternate dictionary learning and a regularized map update
    % The outer loop updates maps and dictionaries. Inside it, dictionary
    % learning repeatedly updates D then C while the maps stay fixed.
    % Backtracking is another loop: it retries one map subproblem with more damping.
    for iteration=1:settings.MaxIterations
        % A. Hold x fixed and learn D/C separately for T1, T2 and density.
        % eta decreases with the outer iteration; its scale always uses initial D/C.
        % In one-step/fixed mode eta is recorded but does not stop the inner loop.
        eta=settings.InnerToleranceScale*initial_dictionary_norm/iteration^settings.InnerDecay;
        % iteration starts at 1, corresponding to k+1 in zero-based notation.
        Dnext=D; Cnext=C; % Work on candidates so a later failure can roll back everything.
        inner_steps=zeros(1,3);
        % Record work spent separately from the final inner-step accuracy.
        inner_changes=zeros(1,3);
        inner_ok=true; % All active map dictionaries must meet the selected inner rule.
        for j=1:3
            if settings.Lambda(j)==0, continue; end % This block is absent from the objective.
            X=create_data_matrix_fast(x(:,:,j),settings.PatchSize,R); % Freeze this map's patches.
            % Warm-start from the dictionary and coefficients of the last outer step.
            inner_options=struct('D0',D(:,:,j),'C0',C(:,:,j),'StopRule',inner_rule, ...
                'MaxIterations',inner_limit,'Tolerance',eta, ...
                'LambdaD',settings.LambdaD,'LambdaC',settings.LambdaC);
            % The inner solver alternates an SVD-based D update and soft-thresholding C.
            [Dnext(:,:,j),Cnext(:,:,j),inner]=dictionary_reconstruction(X,settings.Beta(j),inner_options);
            inner_steps(j)=inner.iterations;
            inner_changes(j)=inner.step_norm(end); % LAST D/C change, used by the inner stop.
            inner_ok=inner_ok && (strcmp(inner_rule,'fixed') || inner.tolerance_met);
        end
        info.last_inner_tolerance=eta;
        info.last_inner_step_norm=inner_changes;
        if ~inner_ok
            info.stop_reason='inner_tolerance_not_met';
            break; % Preserve the previous x, D and C when the inner cap is insufficient.
        end
        % Outer stopping uses the TOTAL change across the entire inner solve.
        dictionary_change=hypot(norm(Dnext(:)-D(:)),norm(Cnext(:)-C(:)));

        % B. Hold candidate D/C fixed and form the quadratic model for x.
        update_patch_targets(Dnext,Cnext);
        [patch_value,patch_gradient]=patch_model(x); % Patch-fit energy and derivative.
        % Compare map candidates against the objective with the SAME new D/C.
        base_value=fdata+.5*settings.Alpha*(x(:)'*(K*x(:)))+patch_value+sparsity_value;
        % The L1 penalty is constant while D/C are fixed: it contributes to
        % this objective but has no derivative with respect to the map variables.
        g=gdata+settings.Alpha*(K*x(:))+patch_gradient; % Data + spatial + patch gradient.
        base_H=Happrox+regularizer_H; % Approximate data curvature + exact prior curvature.
        step_lower=lower(free)-x(free); % Bounds constrain x+d, hence subtract x.
        step_upper=upper(free)-x(free);
        qp_start=zeros(number_of_free,1); % The zero increment is feasible.
        success=false;
        qp_opt=NaN; % Remains NaN when all variables are fixed and no QP is needed.
        for trial=1:settings.MaxBacktracking
            % Keep x, its Bloch derivatives and candidate D/C fixed during
            % retries; only damping and the feasible QP starting step can change.
            step=zeros(3*P,1); % Fixed components always have zero increment.
            if number_of_free>0
                total_trials=total_trials+1;
                if use_exact
                    penalty=regularizer_H+damping*G; % Used by the exact Hessian callback only.
                end
                H=base_H+damping*G; % In exact mode, H supplies only the preconditioner.
                % Minimize g'*d + 0.5*d'*H*d subject to the increment bounds.
                % Empty arguments mean there are no additional linear constraints.
                [free_step,~,qp_exitflag,qp_output]=quadprog(H(free,free),g(free),[],[],[],[], ...
                    step_lower,step_upper,qp_start,qpopts);
                qp_opt=qp_output.firstorderopt;
                info.last_qp_exitflag=qp_exitflag;
                info.last_qp_firstorderopt=qp_opt;
                valid=numel(free_step)==number_of_free && all(isfinite(free_step));
                if valid
                    qp_start=min(step_upper,max(step_lower,free_step)); % Warm-start retries.
                end
                % A positive exit flag alone may mean stagnation; optionally require optimality.
                qp_accurate=~settings.QPRequireOptimality || qp_opt<=settings.QPOptimalityTolerance;
                if ~(valid && qp_exitflag>0 && qp_accurate)
                    damping=damping*4;
                    continue; % Retry the QP before evaluating the nonlinear objective.
                end
                step(free)=free_step; % Embed free increments; fixed components remain zero.
            end
            candidate=min(upper,max(lower,x(:)+step)); % Remove bound roundoff.
            step=candidate-x(:);
            step_squared=max(0,step'*(G*step)); % Squared H01/H1 norm; guard roundoff.
            trial_x=reshape(candidate,n,m,3); % Return to image shape for Bloch evaluation.
            trial_data=data_model(trial_x); % No derivatives needed for a rejected candidate.
            trial_patch=patch_model(trial_x);
            trial_value=trial_data+.5*settings.Alpha*(candidate'*(K*candidate))+trial_patch+sparsity_value;
            % Require nonlinear decrease >= damping*||step||_G^2/4 (up to roundoff).
            % This is the paper's sufficient-decrease test with sigma_3=0.5.
            if isfinite(trial_value) && trial_value<=base_value-.25*damping*step_squared ...
                    +10*eps*max(1,abs(base_value))
                success=true;
                break; % Accept this map candidate and stop increasing damping.
            end
            damping=damping*4; % Failed decrease: penalize large increments more strongly.
        end
        if ~success
            info.stop_reason='backtracking_failed';
            break; % Neither the candidate maps nor the candidate dictionaries are committed.
        end

        % C. Accept maps AND dictionaries together, after every check has passed.
        % Until here, the returned iterate would still be the previous complete triple.
        x=trial_x;
        D=Dnext;
        C=Cnext;
        completed=completed+1;
        info.step_norm(completed)=sqrt(step_squared); % Absolute G norm, as in the paper.
        info.dictionary_step_norm(completed)=dictionary_change; % Across the entire inner solve.
        info.damping(completed)=damping;
        info.backtracking(completed)=trial-1;
        info.qp_firstorderopt(completed)=qp_opt;
        info.inner_tolerance(completed)=eta;
        info.inner_iterations(completed,:)=inner_steps;
        info.inner_step_norm(completed,:)=inner_changes; % LAST inner changes, per map.
        update_data_model(); % Refresh derivatives at the accepted maps for the next QP.
        record_state();      % Evaluation/history only; it does not change x, D or C.
        if settings.Verbose
            fprintf('%3d  J=%.6g  data=%.6g  step=%.2g  DC=%.2g  inner=[%d %d %d]  damping=%.2g\n', ...
                completed,info.objective(completed+1),fdata,sqrt(step_squared), ...
                dictionary_change,inner_steps,damping);
        end
        % D. Stop only when BOTH map change and total dictionary/code change are small.
        % These absolute step tests do not certify joint stationarity.
        if sqrt(step_squared)<=settings.MapTolerance && dictionary_change<=settings.DictionaryTolerance
            info.stop_reason='small_joint_step';
            break;
        end
        damping=max(settings.MinDamping,damping/4); % Try a less damped map step next time.
    end

    %% 7. Return physical maps and trim histories to completed updates
    for name=[state_fields {'background_rms','monitor'}]
        info.(name{1})=info.(name{1})(1:completed+1,:);
    end
    for name=[step_fields {'inner_iterations','inner_step_norm'}]
        info.(name{1})=info.(name{1})(1:completed,:);
    end
    u=offset+x.*scale; % Return physical maps; no background masking here.
    info.iterations=completed;
    info.total_qp_trials=total_trials;
    info.elapsed_seconds=toc(timer);

    %% Helper: data-fit value, gradient and Gauss-Newton curvature at current x
    function update_data_model()
        % Four outputs build the local approximation. Five also cache exact
        % sampled GN products; the local approximation then preconditions them.
        if use_exact
            [fdata,~,gdata,Happrox,data_product]=data_model(x);
        else
            [fdata,~,gdata,Happrox]=data_model(x);
        end
        % Stop before using invalid Bloch values or derivatives in another QP.
        assert(isfinite(fdata) && all(isfinite(gdata)),'Nonfinite data model or gradient.');
    end

    %% Helper: exact Hessian action requested by quadprog (no dense data matrix)
    function products=exact_hessian(~,directions)
        % quadprog supplies free-variable columns. Embed them before the FFT
        % calculation, add prior/damping curvature, then restrict back to free entries.
        full_directions=zeros(3*P,size(directions,2));
        full_directions(free,:)=full(directions); % Embed free increments with zeros in fixed entries.
        products=data_product(full_directions)+penalty*full_directions;
        products=products(free,:); % Restrict the resulting product back to the QP unknowns.
    end

    %% Helper: quantities constant while the map backtracking loop runs
    function update_patch_targets(dict,coeff)
        % targets has one synthesized patch per column. Its values are held
        % constant during a map solve, so D*C need not be recomputed on retries.
        sparsity_value=0;
        for channel=1:3
            if settings.Lambda(channel)==0, continue; end
            targets(:,:,channel)=dict(:,:,channel)*coeff(:,:,channel);
            coefficients=coeff(:,:,channel);
            sparsity_value=sparsity_value+settings.Lambda(channel)*settings.Beta(channel)*sum(abs(coefficients(:)));
        end
    end

    %% Helper: patch-fit energy and its derivative with D/C fixed
    function [fit,gradient]=patch_model(z)
        fit=0;
        % Rejected map candidates need only energy; avoid their gradient work.
        if nargout>=2, gradient=zeros(P,3); end
        for channel=1:3
            if settings.Lambda(channel)==0, continue; end
            Xpatch=create_data_matrix_fast(z(:,:,channel),settings.PatchSize,R);
            residual=Xpatch-targets(:,:,channel);
            % A pixel may occur in many patches; R' adds every corresponding
            % residual contribution when differentiating the patch-fit energy.
            fit=fit+.5*settings.Lambda(channel)*sum(residual(:).^2);
            if nargout>=2
                gradient(:,channel)=settings.Lambda(channel)*(R'*residual(:)); % Adjoint sum, not overlap average.
            end
        end
        if nargout>=2, gradient=gradient(:); end % Same T1/T2/density ordering as x(:).
    end

    %% Helper: record objective terms, stationarity and optional diagnostics
    function record_state()
        k=completed+1; % Index 1 is the initial state; later indices follow accepted steps.
        [fit,patch_gradient]=patch_model(x);
        spatial=.5*settings.Alpha*(x(:)'*(K*x(:)));
        gradient=gdata+settings.Alpha*(K*x(:))+patch_gradient;
        pg=x(:)-min(upper,max(lower,x(:)-gradient)); % Map stationarity with D,C fixed.
        info.objective(k)=fdata+spatial+fit+sparsity_value;
        info.data(k)=fdata;
        info.spatial(k)=spatial;
        info.patch(k)=fit;
        info.sparsity(k)=sparsity_value;
        info.projected_gradient(k)=norm(pg,inf); % This alone does not test D/C stationarity.
        if ~isempty(settings.BodyMask) && any(~settings.BodyMask(:))
            maps=reshape(x,P,3).*settings.Scales+[settings.ReferenceTimes 0];
            info.background_rms(k,:)=sqrt(mean(maps(~settings.BodyMask(:),:).^2,1));
        end
        if ~isempty(settings.MonitorFcn)
            metrics=settings.MonitorFcn(offset+x.*scale); % Evaluation only; never enters the objective.
            if k==1, info.monitor=zeros(settings.MaxIterations+1,numel(metrics)); end
            info.monitor(k,:)=metrics(:).';
        end
        info.time(k)=toc(timer);
    end
end

%% Settings and input checks (not part of the alternating algorithm)
function settings=read_options(options,data,RF,TR,r,zeta,u0)
    % assert(condition,message) stops with that message when condition is false.
    defaults=struct('Alpha',1,'Lambda',[45 45 45],'Beta',[.0045 .0045 .0045], ...
        'PatchSize',[8 8],'Scales',[5000 500 100],'ReferenceTimes',[1000 100], ...
        'H',1,'StepNorm','H1','GaussNewton','approximate', ...
        'LowerBounds',[100 10 0],'UpperBounds',[5500 600 110],'BodyMask',[], ...
        'MaxIterations',40,'MapTolerance',1e-6,'DictionaryTolerance',1e-3, ...
        'InnerStopRule','tolerance','InnerMaxIterations',500, ...
        'InnerToleranceScale',1,'InnerDecay',.75,'LambdaD',1,'LambdaC',1, ...
        'D0',[],'C0',[],'InitialDamping',1,'MinDamping',1e-8,'MaxBacktracking',20, ...
        'QPOptimalityTolerance',1e-8,'QPFunctionTolerance',0, ...
        'QPPCGTolerance',1e-6,'QPRequireOptimality',true,'QPMaxIterations',200, ...
        'Verbose',true,'MonitorFcn',[]);
    % Alpha weights spatial gradients, Lambda weights patch fitting, and
    % Lambda.*Beta weights coefficient sparsity. LambdaD/LambdaC instead
    % penalize changes inside the dictionary-learning subproblems.
    % MapTolerance/DictionaryTolerance govern the OUTER change test;
    % InnerToleranceScale/InnerDecay govern the INNER dictionary accuracy.
    names=fieldnames(options);
    for j=1:numel(names)
        assert(isfield(defaults,names{j}),'Unknown option: %s.',names{j});
        defaults.(names{j})=options.(names{j});
    end
    settings=defaults;
    % Validate the three independent choices: data curvature, step norm,
    % and how much dictionary work is performed per outer update.
    settings.GaussNewton=validatestring(settings.GaussNewton,{'approximate','exact'});
    settings.StepNorm=validatestring(settings.StepNorm,{'H01','H1'});
    settings.InnerStopRule=validatestring(settings.InnerStopRule,{'tolerance','fixed','one-step'});
    [n,m,channels]=size(u0); L=size(data,3);
    % assert stops with its message if a condition fails; these checks keep
    % invalid dimensions, sampling patterns or physical bounds out of the solvers.
    assert(channels==3 && size(data,1)==n && size(data,2)==m, ...
        'Use N-by-M-by-3 physical maps and N-by-M-by-L data.');
    assert(numel(RF)==L && numel(TR)==L && numel(zeta)==L && all(TR(:)>0), ...
        'Supply RF, positive TR and sampling offsets for every frame.');
    assert(r>=1 && r==fix(r) && mod(n,r)==0 && ...
        all(zeta(:)>=1 & zeta(:)<=r & zeta(:)==fix(zeta(:))), 'Invalid row sampling.');
    for name={'Scales','LowerBounds','UpperBounds','Lambda','Beta'}
        settings.(name{1})=reshape(settings.(name{1}),1,3);
    end
    settings.ReferenceTimes=reshape(settings.ReferenceTimes,1,2);
    assert(all(settings.Scales>0) && all(settings.UpperBounds>=settings.LowerBounds) && ...
        all(settings.LowerBounds(1:2)>0) && settings.LowerBounds(3)>=0, 'Invalid scales or physical bounds.');
    assert(settings.Alpha>=0 && all(settings.Lambda>=0) && all(settings.Beta>=0) && settings.H>0, ...
        'Regularization weights must be nonnegative and H positive.');
    assert(settings.InitialDamping>0 && settings.MinDamping>0 && settings.InnerToleranceScale>0 && ...
        settings.InnerDecay>0 && settings.LambdaD>0 && settings.LambdaC>0, 'Use positive damping and inner settings.');
    assert(settings.MaxIterations>=1 && settings.MaxIterations==fix(settings.MaxIterations) && ...
        settings.InnerMaxIterations>=1 && settings.InnerMaxIterations==fix(settings.InnerMaxIterations) && ...
        settings.MaxBacktracking>=1 && settings.MaxBacktracking==fix(settings.MaxBacktracking), ...
        'Iteration limits must be positive integers.');
    assert(isempty(settings.BodyMask) || isequal(size(settings.BodyMask),[n m]), 'BodyMask must match the grid.');
end

%% Allocate output histories once; these arrays do not drive the updates
function [info,state_fields,step_fields]=initialize_history(max_iterations)
    % The extra row in state histories stores the starting point. The three
    % columns of inner histories correspond to T1, T2 and density dictionaries.
    state_fields={'objective','data','spatial','patch','sparsity','projected_gradient','time'};
    step_fields={'step_norm','dictionary_step_norm','damping','backtracking', ...
        'qp_firstorderopt','inner_tolerance'};
    for name=state_fields, info.(name{1})=zeros(max_iterations+1,1); end
    for name=step_fields, info.(name{1})=zeros(max_iterations,1); end
    info.inner_iterations=zeros(max_iterations,3);
    info.inner_step_norm=zeros(max_iterations,3);
    info.background_rms=NaN(max_iterations+1,3);
    info.monitor=zeros(max_iterations+1,0); % No additional histories unless requested.
    info.stop_reason='max_iterations'; info.last_qp_exitflag=NaN;
    info.last_qp_firstorderopt=NaN;
end
