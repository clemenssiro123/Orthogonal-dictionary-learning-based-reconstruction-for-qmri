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
%   The objective, with theta=(u-[ReferenceTimes,0])./Scales, is
%
%     h^2/2*||Q(u)-data||^2 + Alpha/2*theta'*K_grad*theta
%       + sum_j Lambda(j)*(0.5*||R*theta_j-D_j*C_j||_F^2 + Beta(j)*||C_j||_1).
%
%   R extracts overlapping patches; each column of C represents one patch.
%   Each D(:,:,j) is orthogonal. K_grad=h^2*diag(-LAP,-LAP,-LAP) penalizes
%   spatial gradients, with zero exterior values for theta. Beta is INSIDE Lambda*(...).
%   This is the discrete P2/Appendix A convention used by this project.
%
%   Algorithm 2 notation (corrected_thesis_version.pdf, p. 18):
%     k starts at zero; theta_k represents the paper's u_k in scaled coordinates.
%     z_k=(D_k,C_k); the candidate z_{k+1} is (D_kp1,C_kp1).
%     eta_k is the inner tolerance; epsilon_1/epsilon_2 are outer tolerances.
%     Lambda_0 is the initial trial parameter; lambda_min is underline{lambda}.
%     Lambda_k starts backtracking; lambda_trial is its running lambda.
%     lambda_k is the accepted parameter; tau=4 and sigma_BT=0.5.
%     theta_lambda represents u(lambda); delta_theta=theta_lambda(:)-theta_k(:).
%   settings.Lambda contains PATCH weights, distinct from Lambda_k and lambda_k.
%   G_step=K_grad gives the discrete H01 norm in Algorithm 2; the H1 option
%   adds h^2*I. Exact GN uses the sampled linearized data model; approximate GN
%   uses the paper's local Gram approximation. QP tolerances, finite iteration
%   caps, bound clipping and a roundoff allowance are numerical safeguards.
%
%   InnerStopRule: 'tolerance' = nested; 'one-step' = one D/C pair per map;
%   'fixed' = a prescribed number of pairs. The map update is the same.
%   The outer stop requires BOTH small map and small total D/C changes.
%   A failed inner solve or map update returns the last COMPLETE iterate.
%   BodyMask/MonitorFcn are diagnostics only; all pixels remain variables.
%   Each LowerBounds entry must be strictly below its UpperBounds entry.
%   MATLAB R2021b + Optimization Toolbox. See demo_compare_all_reconstructions.
%
%   Array sizes: data is N-by-M-by-L; u/u0 are N-by-M-by-3 physical maps.
%   For A=prod(PatchSize) and Q overlapping patches, D is A-by-A-by-3 and
%   C is A-by-Q-by-3. Channel order is always T1,T2,rho. Returned D/C describe
%   patches of scaled deviations theta_k, not physical u. State histories include
%   the initial point; step histories describe only completed outer updates.

    if nargin<7, options=struct(); end
    root=fileparts(mfilename('fullpath'));
    addpath(fullfile(root,'helper_functions'),fullfile(root,'bloch_simulations'));

    %% 1. Settings: choose nested or one-step dictionary learning
    % Fill omitted options and check the input geometry before allocating the
    % reconstruction arrays. All three dictionary-update modes use the same
    % map QP and acceptance test; only the work inside the z-step changes.
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
    P=n*m; % Pixels per parameter map; the map QP has 3*P unknowns.
    RF=RF(:).'; TR=TR(:).'; zeta=zeta(:).'; % .' preserves complex RF phases.

    %% 2. Scale physical maps and translate the physical bounds
    % The optimization variables are dimensionless deviations from the
    % reference maps. Data evaluation converts them back to physical units;
    % the spatial and dictionary penalties act directly on these deviations.
    timer=tic;
    scale=reshape(settings.Scales,1,1,3); % One scale per channel, broadcast over pixels.
    offset=reshape([settings.ReferenceTimes 0],1,1,3); % Origin [T1ref,T2ref,0].
    theta_k=(u0-offset)./scale; % Dimensionless deviations; these may be negative.
    % Bounds must be expressed in the same coordinates as the QP unknowns.
    lower=repmat((reshape(settings.LowerBounds,1,1,3)-offset)./scale,n,m,1);
    upper=repmat((reshape(settings.UpperBounds,1,1,3)-offset)./scale,n,m,1);
    lower=lower(:); upper=upper(:); % Map-major ordering, identical to theta_k(:).
    % Stop here if the starting maps are nonfinite or physically infeasible.
    assert(all(isfinite(theta_k(:)) & theta_k(:)>=lower & theta_k(:)<=upper),'Initial maps must satisfy bounds.');

    %% 3. Build patch extraction, spatial regularization and the step metric
    % R acts on one vectorized map; counts gives each pixel's patch coverage.
    % Its geometry depends only on image size, patch size and stride, so R
    % is built once here and reused for all maps and all outer/inner trials.
    % Later calls to create_data_matrix_fast extract new values with this R.
    [R,number_of_patches,counts]=create_R_matrix([n m],settings.PatchSize,1);
    patch_area=prod(settings.PatchSize); % Number of entries in a single patch.
    % theta_k(:)'*K_grad*theta_k(:) discretizes the summed gradient energy.
    % kron(speye(3),...) places one zero-Dirichlet Laplacian on each map's
    % diagonal block. The h^2 factor supplies the 2-D spatial integration
    % weight; the Laplacian itself already contains the factor 1/h^2.
    K_grad=settings.H^2*kron(speye(3),-dirichlet_laplacian_matrix(n,m,settings.H));
    G_step=K_grad; % Algorithm 2: ||delta_theta||_H01^2 = delta_theta'*G_step*delta_theta.
    % The optional full H1 norm also measures the increment values themselves.
    if strcmp(settings.StepNorm,'H1'), G_step=G_step+settings.H^2*speye(3*P); end
    patch_diagonal=counts*settings.Lambda; % One overlap-count vector per parameter.
    % Outer product: entry (p,j) is counts(p)*Lambda(j), a P-by-3 array.
    % R'*R is diagonal, so the patch-fit Hessian needs only overlap counts.
    regularizer_H=settings.Alpha*K_grad+spdiags(patch_diagonal(:),0,3*P,3*P);

    %% 4. Initialize an orthogonal dictionary and sparse coefficients per map
    % Supplied D0/C0 are used as starting values. Missing values are built
    % from the initial scaled maps: an uncentered patch second-moment matrix
    % gives the orthogonal basis, and soft-thresholding gives its sparse codes.
    D_k=settings.D0;
    C_k=settings.C0;
    if isempty(D_k), D_k=zeros(patch_area,patch_area,3); end
    if isempty(C_k), C_k=zeros(patch_area,number_of_patches,3); end
    assert(isequal(size(D_k),[patch_area patch_area 3]) && ...
        isequal(size(C_k),[patch_area number_of_patches 3]), 'D0 or C0 has incompatible size.');
    for j=1:3
        X=create_data_matrix_fast(theta_k(:,:,j),settings.PatchSize,R); % One patch per column.
        % PCA of the starting patches supplies an orthogonal basis if D0 is absent.
        if isempty(settings.D0), [D_k(:,:,j),~,~]=svd(X*X','econ'); end
        % Orthogonality is required for the closed-form D/C update formulas.
        assert(norm(D_k(:,:,j)'*D_k(:,:,j)-eye(patch_area),'fro')<1e-8*sqrt(patch_area), ...
            'Each initial dictionary must be orthogonal.');
        % Project onto the basis and soft-threshold to initialize sparse codes.
        if isempty(settings.C0), C_k(:,:,j)=prox_c(settings.Beta(j),D_k(:,:,j)'*X,'soft'); end
    end
    z0_norm=hypot(norm(D_k(:)),norm(C_k(:))); % Fixed scale in eta_k.
    % ||z_0||_Z stays fixed even when the norms of later D/C iterates change.

    %% 5. Prepare data derivatives, the bounded QP solver and output histories
    % This function handle binds the fixed measurements, pulse sequence and
    % coordinate scaling. A call with one output evaluates only the data
    % objective; requesting more outputs also builds the needed derivatives.
    data_model=@(z) qmri_data_objective(z,data,RF,TR,r,zeta,settings.Scales,settings.H,settings.ReferenceTimes);
    penalty=[]; % Prior/damping matrix used by the exact QP Hessian callback.
    [targets,sparsity_value]=patch_targets(D_k,C_k);
    [fdata,gdata,Happrox,data_product]=evaluate_data_model(theta_k);
    qpopts=optimoptions('quadprog','Algorithm','trust-region-reflective','Display','off', ...
        'OptimalityTolerance',settings.QPOptimalityTolerance,'FunctionTolerance',settings.QPFunctionTolerance, ...
        'TolPCG',settings.QPPCGTolerance,'PrecondBandWidth',inf,'MaxIterations',settings.QPMaxIterations); 
    % Exact mode uses FFT products for the data Hessian; the sparse matrix preconditions it.
    if use_exact, qpopts.HessianMultiplyFcn=@exact_hessian; end

    [info,state_fields,step_fields]=initialize_history(settings.MaxIterations);
    completed=0;  % Count only complete, accepted outer updates.
    total_trials=0; % Count every map QP attempt, including rejected attempts.
    % Algorithm 2, line 1. Keep the public option/history names compatible.
    epsilon_1=settings.MapTolerance;
    epsilon_2=settings.DictionaryTolerance;
    tau=4; sigma_BT=0.5; % Backtracking multiplier and sufficient-decrease fraction.
    lambda_min=settings.MinDamping; % Paper: underline{lambda}.
    Lambda_0=max(settings.InitialDamping,lambda_min); % Enforce Lambda_0 >= lambda_min.
    lambda_k=NaN; % No accepted parameter yet; read only after a completed update.
    record_state(); % History entry 1 describes the supplied starting maps.

    %% 6. Algorithm 2: alternate the z-step and the regularized u-step
    % The outer loop updates maps and dictionaries. Inside it, dictionary
    % learning repeatedly updates D then C while the maps stay fixed.
    % Backtracking retries one map subproblem with a larger lambda_trial.
    for k=0:(settings.MaxIterations-1)
        % A. Algorithm 2, lines 4-5: z_k=(D_k,C_k) -> z_{k+1}=(D_kp1,C_kp1).
        % Hold theta_k fixed and learn D/C separately for T1, T2 and density.
        % eta_k decreases with k; its scale always uses the initial z_0.
        % In one-step/fixed mode eta_k is recorded but does not stop the inner loop.
        eta_k=settings.InnerToleranceScale*z0_norm/(k+1)^settings.InnerDecay;
        D_kp1=D_k; C_kp1=C_k; % Work on candidates so a later failure can roll back everything.
        inner_steps=zeros(1,3);
        % Record work spent separately from the final inner-step accuracy.
        inner_changes=zeros(1,3);
        inner_ok=true; % All active map dictionaries must meet the selected inner rule.
        for j=1:3
            if settings.Lambda(j)==0, continue; end % This block is absent from the objective.
            X=create_data_matrix_fast(theta_k(:,:,j),settings.PatchSize,R); % Freeze this map's patches.
            % Warm-start the inner solve from z_k, the last accepted dictionary/code pair.
            inner_options=struct('D0',D_k(:,:,j),'C0',C_k(:,:,j),'StopRule',inner_rule, ...
                'MaxIterations',inner_limit,'Tolerance',eta_k, ...
                'LambdaD',settings.LambdaD,'LambdaC',settings.LambdaC);
            % The inner solver alternates an SVD-based D update and soft-thresholding C.
            [D_kp1(:,:,j),C_kp1(:,:,j),inner]=dictionary_reconstruction(X,settings.Beta(j),inner_options);
            inner_steps(j)=inner.iterations;
            inner_changes(j)=inner.step_norm(end); % LAST inner D/C change, per map.
            inner_ok=inner_ok && (strcmp(inner_rule,'fixed') || inner.tolerance_met);
        end
        info.last_inner_tolerance=eta_k;
        info.last_inner_step_norm=inner_changes;
        if ~inner_ok
            info.stop_reason='inner_tolerance_not_met';
            break; % Preserve the previous theta_k, D_k and C_k when the inner cap is insufficient.
        end
        % ||z_{k+1}-z_k||_Z: TOTAL change across the entire inner solve.
        z_step_norm=hypot(norm(D_kp1(:)-D_k(:)),norm(C_kp1(:)-C_k(:)));

        % B. Algorithm 2, lines 6-11: u-step with candidate z_{k+1} fixed.
        % Equation (38): choose Lambda_k at the START of the current u-step.
        % Before this iteration accepts a step, lambda_k holds lambda_{k-1}.
        % Try one damping level below the preceding accepted value, subject
        % to the lower safeguard. If this trial fails, multiplication by tau
        % returns to the preceding value when the safeguard was inactive.
        if k==0
            Lambda_k=Lambda_0;
        else
            Lambda_k=max(lambda_min,lambda_k/tau);
        end
        lambda_trial=Lambda_k; % Paper's running lambda; increase it after rejection.
        % Refresh the shared patch targets for z_{k+1}. Both the baseline
        % objective and every trial objective must use this same dictionary
        % block; the data linearization still belongs to the current theta_k.
        [targets,sparsity_value]=patch_targets(D_kp1,C_kp1);
        [patch_value,patch_gradient]=patch_model(theta_k); % Patch-fit energy and derivative.
        % J_k_zkp1 = J(u_k,z_{k+1}), using the SAME new D/C as every candidate.
        J_k_zkp1=fdata+.5*settings.Alpha*(theta_k(:)'*(K_grad*theta_k(:)))+patch_value+sparsity_value;
        % The L1 penalty is constant while D_kp1/C_kp1 are fixed: it contributes to
        % this objective but has no derivative with respect to the map variables.
        g=gdata+settings.Alpha*(K_grad*theta_k(:))+patch_gradient; % Data + spatial + patch gradient.
        base_H=Happrox+regularizer_H; % Approximate data curvature + exact prior curvature.
        step_lower=lower-theta_k(:); % Bounds constrain theta_k+d, hence subtract theta_k.
        step_upper=upper-theta_k(:);
        % The QP unknown is an increment, so starting at zero means starting
        % at the current feasible maps. This solver guess is separate from
        % Lambda_k, which controls the quadratic penalty on that increment.
        qp_start=zeros(3*P,1); % The zero increment is feasible.
        success=false;
        for trial=1:settings.MaxBacktracking
            % Keep theta_k, its Bloch derivatives and candidate D_kp1/C_kp1 fixed.
            % Only lambda_trial and the feasible QP starting increment can change.
            total_trials=total_trials+1;
            if use_exact
                penalty=regularizer_H+lambda_trial*G_step; % Used by the exact Hessian callback only.
            end
            H=base_H+lambda_trial*G_step; % In exact mode, H supplies only the preconditioner.
            % Equation (39), written in delta_theta = theta_lambda - theta_k.
            % Exact mode gets curvature from exact_hessian; H preconditions it.
            [delta_theta,~,qp_exitflag,qp_output]=quadprog(H,g,[],[],[],[], ...
                step_lower,step_upper,qp_start,qpopts);
            qp_opt=qp_output.firstorderopt;
            info.last_qp_exitflag=qp_exitflag;
            info.last_qp_firstorderopt=qp_opt;
            valid=numel(delta_theta)==3*P && all(isfinite(delta_theta));
            if valid
                qp_start=min(step_upper,max(step_lower,delta_theta)); % Warm-start retries.
            end
            % A positive exit flag alone may mean stagnation; optionally require optimality.
            % These checks concern the quadratic subproblem only. A valid,
            % accurate QP solution must still pass the nonlinear decrease
            % test below before it can replace the current reconstruction.
            qp_accurate=~settings.QPRequireOptimality || qp_opt<=settings.QPOptimalityTolerance;
            if ~(valid && qp_exitflag>0 && qp_accurate)
                lambda_trial=tau*lambda_trial; % Additional numerical QP safeguard.
                continue; % Retry the QP before evaluating the nonlinear objective.
            end
            candidate=min(upper,max(lower,theta_k(:)+delta_theta)); % Remove bound roundoff.
            delta_theta=candidate-theta_k(:);
            theta_step_squared=max(0,delta_theta'*(G_step*delta_theta)); % Squared H01/H1 norm; guard roundoff
            theta_lambda=reshape(candidate,n,m,3); % Paper's u(lambda), in scaled coordinates.
            trial_data=data_model(theta_lambda); % No derivatives needed for a rejected candidate.
            trial_patch=patch_model(theta_lambda);
            J_lambda_zkp1=trial_data+.5*settings.Alpha*(candidate'*(K_grad*candidate))+trial_patch+sparsity_value;
            % Algorithm 2, line 8 / equation (40):
            % J(u(lambda),z_{k+1}) <= J(u_k,z_{k+1}) - sigma_BT*lambda/2*||delta_theta||_G^2.
            % Evaluate the original objective, including its constant L1
            % term. Damping belongs to the QP model and the required decrease;
            % it is not an extra term in J_lambda_zkp1. The eps allowance
            % prevents rejection caused solely by objective roundoff.
            if isfinite(J_lambda_zkp1) && J_lambda_zkp1<=J_k_zkp1-(sigma_BT/2)*lambda_trial*theta_step_squared ...
                    +10*eps*max(1,abs(J_k_zkp1))
                success=true;
                break; % Accept this map candidate and stop increasing lambda_trial.
            end
            lambda_trial=tau*lambda_trial; % Algorithm 2, line 9: lambda <- tau*lambda.
        end
        if ~success
            info.stop_reason='backtracking_failed';
            break; % Neither the candidate maps nor the candidate dictionaries are committed.
        end

        % C. Algorithm 2, line 12: accept lambda_k and u_{k+1}=u(lambda_k).
        % Commit maps AND dictionaries together after every check has passed.
        % Until here, the returned iterate would still be the previous complete triple.
        lambda_k=lambda_trial;
        theta_k=theta_lambda;
        D_k=D_kp1;
        C_k=C_kp1;
        completed=completed+1;
        info.step_norm(completed)=sqrt(theta_step_squared); % Absolute G_step norm, as in the paper.
        info.dictionary_step_norm(completed)=z_step_norm; % Across the entire inner solve.
        info.damping(completed)=lambda_k; % Existing history field stores accepted lambda_k.
        info.backtracking(completed)=trial-1;
        info.qp_firstorderopt(completed)=qp_opt;
        info.inner_tolerance(completed)=eta_k;
        info.inner_iterations(completed,:)=inner_steps;
        info.inner_step_norm(completed,:)=inner_changes; % LAST inner changes, per map.
        [fdata,gdata,Happrox,data_product]=evaluate_data_model(theta_k);
        record_state();      % Evaluation/history only; it does not change theta_k, D_k or C_k.
        if settings.Verbose
            fprintf('%3d  J=%.6g  data=%.6g  step=%.2g  DC=%.2g  inner=[%d %d %d]  lambda_k=%.2g\n', ...
                completed,info.objective(completed+1),fdata,sqrt(theta_step_squared), ...
                z_step_norm,inner_steps,lambda_k);
        end
        % D. Algorithm 2, line 14: both changes must meet epsilon_1/epsilon_2.
        % Use <= at the numerical tolerance boundary; these tests do not certify stationarity.
        if sqrt(theta_step_squared)<=epsilon_1 && z_step_norm<=epsilon_2
            info.stop_reason='small_joint_step';
            break;
        end
        % The for loop advances k (Algorithm 2, line 13); the next u-step forms Lambda_k.
    end

    %% 7. Return physical maps and trim histories to completed updates
    % Failed attempts can update the last_* diagnostics, but do not add an
    % accepted state. Keep one state row for the initial point plus one per
    % completed update; step histories contain only the completed updates.
    for name=[state_fields {'background_rms','monitor'}]
        info.(name{1})=info.(name{1})(1:completed+1,:);
    end
    for name=[step_fields {'inner_iterations','inner_step_norm'}]
        info.(name{1})=info.(name{1})(1:completed,:);
    end
    u=offset+theta_k.*scale; % Return physical maps; no background masking here.
    D=D_k; C=C_k; % Preserve the public output names and the last complete iterate.
    info.iterations=completed;
    info.total_qp_trials=total_trials;
    info.elapsed_seconds=toc(timer);

    % The five helpers below are NESTED functions: they can read variables
    % from this reconstruction call without receiving them as arguments.
    % In particular, exact_hessian reads data_product/penalty, patch_model
    % reads targets, and record_state writes info. Their shared quantities
    % are refreshed explicitly by the outer algorithm at the places above.
    % The helper argument z denotes scaled maps, not the paper's z=(D,C).

    %% Helper: data-fit value, gradient and Gauss-Newton curvature at current theta_k
    function [value,gradient,local_hessian,product]=evaluate_data_model(z)
        % Input z is n-by-m-by-3 in the same scaled coordinates as theta_k.
        % This wrapper requests the data quantities needed for the next QP:
        %   value         scalar h^2/2*||Q(physical(z))-data||^2 on measured rows;
        %   gradient      real 3*P column vector, in T1/T2/density order;
        %   local_hessian sparse 3*P-by-3*P local Gram approximation;
        %   product       handle applying the exact sampled GN matrix, or [].
        % These outputs contain data terms only. The caller adds the spatial
        % and patch derivatives, and adds damping separately for each trial.
        % Rebuild derivatives only at the initial or accepted maps. Backtracking
        % reuses this linearization and calls data_model directly for trial values.
        % In exact mode the fifth output captures the Bloch sensitivities at z.
        % The ignored second output is a structure of objective components.
        product=[];
        if use_exact
            [value,~,gradient,local_hessian,product]=data_model(z);
        else
            [value,~,gradient,local_hessian]=data_model(z);
        end
        assert(isfinite(value) && all(isfinite(gradient)),'Nonfinite data model or gradient.');
    end

    %% Helper: exact Hessian action requested by quadprog (no dense data matrix)
    function products=exact_hessian(~,directions)
        % quadprog can request several Hessian-vector products at once;
        % each column of directions is one 3*P-dimensional map increment.
        % The ignored first argument is the supplied H matrix, whose role
        % here is preconditioning. The actual curvature is applied below:
        %   H_exact*V = h^2*real(J_data'*J_data*V) + penalty*V.
        % data_product holds the linearization at theta_k throughout retries;
        % penalty=regularizer_H+lambda_trial*G_step changes with each trial.
        % Thus this callback needs no new Bloch simulation. Only the direction
        % columns are made full for FFT operations; no dense Hessian is formed.
        directions=full(directions); % FFT products reshape each column into three maps.
        products=data_product(directions)+penalty*directions;
    end

    %% Helper: quantities constant while the map backtracking loop runs
    function [patches,sparsity]=patch_targets(dict,coeff)
        % Input dict contains three patch_area-by-patch_area dictionaries;
        % coeff contains three patch_area-by-number_of_patches code matrices.
        % Output patches has the same size as coeff, with one reconstructed
        % patch per column: patches(:,:,j)=dict(:,:,j)*coeff(:,:,j).
        % sparsity is the scalar sum_j Lambda(j)*Beta(j)*||coeff(:,:,j)||_1,
        % where the L1 norm sums absolute values of ALL coefficient entries.
        % Both outputs stay fixed while only the maps change. The caller
        % stores them in targets/sparsity_value before evaluating patch_model.
        % Initial history uses D_k/C_k; map backtracking uses D_kp1/C_kp1.
        patches=zeros(size(coeff));
        sparsity=0;
        for channel=1:3
            % A zero patch weight removes both the fit and sparsity terms
            % for this channel, so its initialized zero target is sufficient.
            if settings.Lambda(channel)==0, continue; end
            patches(:,:,channel)=dict(:,:,channel)*coeff(:,:,channel);
            coefficients=coeff(:,:,channel);
            sparsity=sparsity+settings.Lambda(channel)*settings.Beta(channel)*sum(abs(coefficients(:)));
        end
    end

    %% Helper: patch-fit energy and its derivative with the current targets fixed
    function [fit,gradient]=patch_model(z)
        % For scaled maps z and the shared targets, compute the smooth term
        %   fit = sum_j Lambda(j)/2 * ||patches(z_j)-targets_j||_F^2.
        % Neither coefficient sparsity nor spatial regularization is included;
        % those terms are added separately when forming the full objective.
        % With two outputs, also return its real 3*P map gradient. With one
        % output, skip the gradient allocation and adjoint multiplication.
        % nargout is the number of outputs requested from this helper call.
        fit=0;
        % Rejected map candidates need only energy; avoid their gradient work.
        if nargout>=2, gradient=zeros(P,3); end
        for channel=1:3
            if settings.Lambda(channel)==0, continue; end
            % Each column is one vectorized patch. Vectorizing the resulting
            % residual keeps the ordering used by the extraction matrix R.
            Xpatch=create_data_matrix_fast(z(:,:,channel),settings.PatchSize,R);
            residual=Xpatch-targets(:,:,channel);
            % A pixel may occur in many patches; R' adds every corresponding
            % residual contribution when differentiating the patch-fit energy.
            % In vector form this is Lambda(j)*R'*(R*z_j(:)-targets_j(:)).
            % Dividing by overlap counts would change this derivative: each
            % occurrence contributes separately to the objective's patch sum.
            fit=fit+.5*settings.Lambda(channel)*sum(residual(:).^2);
            if nargout>=2
                gradient(:,channel)=settings.Lambda(channel)*(R'*residual(:)); % Adjoint sum, not overlap average.
            end
        end
        if nargout>=2, gradient=gradient(:); end % Same T1/T2/density ordering as theta_k(:).
    end

    %% Helper: record objective terms and optional diagnostics
    function record_state()
        % Read the current accepted maps and write one row into shared info.
        % Before calling, the outer code must ensure that fdata belongs to
        % theta_k and targets/sparsity_value belong to its accepted D_k/C_k.
        % This helper evaluates the spatial and patch terms, then combines
        % all four objective contributions; it does not perform an update.
        % completed=0 stores the initial state in row 1. After an acceptance,
        % row completed+1 stores the new state and step row completed describes
        % the transition into it. Rejected trials never call record_state.
        history_index=completed+1; % Index 1 is the initial state; later indices follow accepted steps.
        fit=patch_model(theta_k);
        spatial=.5*settings.Alpha*(theta_k(:)'*(K_grad*theta_k(:)));
        info.objective(history_index)=fdata+spatial+fit+sparsity_value;
        info.data(history_index)=fdata;
        info.spatial(history_index)=spatial;
        info.patch(history_index)=fit;
        info.sparsity(history_index)=sparsity_value;
        % Background RMS is a diagnostic in physical units, computed for
        % pixels outside BodyMask. It does not mask variables in the solver.
        if ~isempty(settings.BodyMask) && any(~settings.BodyMask(:))
            maps=reshape(theta_k,P,3).*settings.Scales+[settings.ReferenceTimes 0];
            info.background_rms(history_index,:)=sqrt(mean(maps(~settings.BodyMask(:),:).^2,1));
        end
        % A user-supplied monitor receives physical maps and returns a fixed
        % number of metrics. Its first result determines the history width;
        % these values are recorded but never used by the acceptance test.
        if ~isempty(settings.MonitorFcn)
            metrics=settings.MonitorFcn(offset+theta_k.*scale); % Evaluation only; never enters the objective.
            if history_index==1, info.monitor=zeros(settings.MaxIterations+1,numel(metrics)); end
            info.monitor(history_index,:)=metrics(:).';
        end
        info.time(history_index)=toc(timer);
    end
end

%% Settings and input checks (not part of the alternating algorithm)
function settings=read_options(options,data,RF,TR,r,zeta,u0)
    % This is a separate local function, so it receives everything it needs
    % through arguments rather than sharing the reconstruction workspace.
    % Start from defaults and replace only fields supplied by the caller.
    % Reject unknown names to expose misspelled options instead of silently
    % running with an unintended default. Return the merged, normalized struct.
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
    % QP tolerances apply to each bounded map subproblem. The iteration limits
    % cap three different loops: outer updates, inner D/C pairs and QP retries.
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
    % Check grid dimensions, sampling and physical bounds before simulation.
    % One RF pulse, repetition time and sampling offset are needed per frame.
    % Requiring r to divide n matches the regular Cartesian row-sampling grid.
    assert(channels==3 && size(data,1)==n && size(data,2)==m, ...
        'Use N-by-M-by-3 physical maps and N-by-M-by-L data.');
    assert(numel(RF)==L && numel(TR)==L && numel(zeta)==L && all(TR(:)>0), ...
        'Supply RF, positive TR and sampling offsets for every frame.');
    assert(r>=1 && r==fix(r) && mod(n,r)==0 && ...
        all(zeta(:)>=1 & zeta(:)<=r & zeta(:)==fix(zeta(:))), 'Invalid row sampling.');
    % Use row vectors consistently: counts*Lambda is an outer product, and
    % physical/scaled coordinate conversion expects one entry per map.
    for name={'Scales','LowerBounds','UpperBounds','Lambda','Beta'}
        settings.(name{1})=reshape(settings.(name{1}),1,3);
    end
    settings.ReferenceTimes=reshape(settings.ReferenceTimes,1,2);
    % Positive scales preserve bound ordering under the coordinate change.
    % Strictly positive relaxation-time bounds keep all feasible Bloch
    % evaluations in their physical domain; density may be zero.
    assert(all(settings.UpperBounds>settings.LowerBounds), ...
        'Each LowerBounds entry must be strictly below its UpperBounds entry.');
    assert(all(settings.Scales>0) && ...
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
    % This separate local function allocates enough rows for the iteration
    % cap. Preallocation avoids repeatedly growing arrays during reconstruction.
    % The returned field-name lists tell the caller which histories to trim:
    %   state_fields: max_iterations+1 rows, including the initial objective;
    %   step_fields:  max_iterations rows, one per accepted outer transition.
    % The three columns of inner histories correspond to T1, T2 and density.
    state_fields={'objective','data','spatial','patch','sparsity','time'};
    step_fields={'step_norm','dictionary_step_norm','damping','backtracking', ...
        'qp_firstorderopt','inner_tolerance'};
    for name=state_fields, info.(name{1})=zeros(max_iterations+1,1); end
    for name=step_fields, info.(name{1})=zeros(max_iterations,1); end
    info.inner_iterations=zeros(max_iterations,3);
    info.inner_step_norm=zeros(max_iterations,3);
    % NaN marks unavailable background diagnostics. Monitor width is unknown
    % until its first evaluation, so it initially has zero columns.
    info.background_rms=NaN(max_iterations+1,3);
    info.monitor=zeros(max_iterations+1,0); % No additional histories unless requested.
    % Exhausting the outer loop leaves this default stop reason in place;
    % an earlier stopping condition replaces it. The last-QP fields remain
    % NaN until a QP is attempted, and can retain a rejected trial's diagnostics.
    info.stop_reason='max_iterations'; info.last_qp_exitflag=NaN;
    info.last_qp_firstorderopt=NaN;
end
