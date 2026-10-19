function [D, C, info] = dictionary_reconstruction(X, beta, options)
%DICTIONARY_RECONSTRUCTION Learn an orthogonal dictionary for fixed patches.
%   [D,C,INFO] = dictionary_reconstruction(X,BETA,OPTIONS) approximately solves
%
%       min  0.5*||X - D*C||_F^2 + beta*sum(abs(C(:))),   D'*D = I.
%
%   X is a real, full, double K-by-M matrix: one image patch per column.
%   D is K-by-K and C is K-by-M. The input patches stay fixed during the
%   solve. No image normalization, mean removal, or patch extraction is
%   performed here. For qMRI, extract patches of the scaled deviations
%   x=(u-[ReferenceTimes,0])./Scales, as in nested_dictionary_reconstruction.
%   Reconstructed patches are D*C; average their overlaps to obtain an image.
%
%   This is the proximal alternating method in Appendix A of
%   corrected_thesis_version.tex, equations dict_learning_PD/PC:
%       Dnew = U*V', where [U,~,V] = svd(X*C' + LambdaD*D);
%       Cnew = soft((Dnew'*X + LambdaC*C)/(1 + LambdaC),
%                   beta/(1 + LambdaC)).
%   The positive proximal weights LambdaD/LambdaC are constant in this call.
%   beta is the coefficient INSIDE the objective above. For an outer term
%   lambda/2*||X-D*C||_F^2 + mu*||C||_1, pass beta = mu/lambda (lambda > 0).
%
%   OPTIONS is an optional scalar struct with these fields:
%     StopRule       'fixed' (default) or 'tolerance'.
%     MaxIterations  Positive integer, default 20. In fixed mode, perform
%                    exactly this many D/C update pairs. In tolerance mode,
%                    this is a hard cap, even if the accuracy is not reached.
%     Tolerance      Absolute iterate-change tolerance eta, default 1e-3.
%                    Tolerance mode stops after an update pair when
%                    sqrt(||Dnew-D||_F^2 + ||Cnew-C||_F^2) <= eta.
%                    This is the paper's stopping rule, not an image-error
%                    or global-optimality guarantee. It depends on patch
%                    count and scaling. For an outer iteration, supply eta_k.
%     LambdaD        Positive dictionary proximal weight, default 1.
%     LambdaC        Positive coefficient proximal weight, default 1.
%     D0             Initial orthogonal K-by-K dictionary. Default: PCA
%                    basis of X, without centering or access to ground truth.
%     C0             Initial K-by-M coefficients. Default: soft(D0'*X,beta).
%                    Supply both D0 and C0 to warm-start a later qMRI step.
%
%   INFO fields:
%     iterations, stop_reason ('fixed_iterations', 'tolerance', or
%     'max_iterations'), tolerance_met, elapsed_seconds;
%     objective and nonzero_fraction: initial value followed by each update;
%     step_norm: one entry per completed update pair (the stopping quantity).
%
%   Example: exactly 10 steps, then continue to an absolute tolerance.
%     [D,C,a] = dictionary_reconstruction(X,0.03,struct('MaxIterations',10));
%     opts = struct('StopRule','tolerance','Tolerance',1e-2, ...
%                   'MaxIterations',100,'D0',D,'C0',C);
%     [D,C,b] = dictionary_reconstruction(X,0.03,opts);
%
%   MATLAB R2021b; no additional toolbox required.

    if nargin < 3
        options = struct();
    end

    %% 1. Read settings and validate the fixed patch data
    options = read_options(X,beta,options);

    %% 2. Initialize an orthogonal basis D and its sparse patch coefficients C
    % K is the number of pixels per patch; M is the number of training patches.
    % There are K atoms, so D is square and its transpose is its inverse.
    timer = tic;
    [K, M] = size(X);
    if isempty(options.D0)
        % Only a K-by-K covariance matrix is needed for PCA initialization.
        [D, ~, ~] = svd(X*X', 'econ');
    else
        D = options.D0;
        validateattributes(D, {'double'}, ...
            {'real','finite','nonsparse','size',[K K]}, mfilename, 'D0');
        if norm(D'*D - eye(K), 'fro') > 1e-8*sqrt(K)
            error('dictionary_reconstruction:NonorthogonalD0', ...
                'D0 must satisfy D0''*D0 = I.');
        end
    end
    projected = D'*X; % Express all fixed patches in the current dictionary basis.
    if isempty(options.C0)
        % With D fixed and orthogonal, the coefficient minimizer is soft(D'*X,beta).
        C = sign(projected).*max(abs(projected) - beta, 0);
    else
        C = options.C0;
        validateattributes(C, {'double'}, ...
            {'real','finite','nonsparse','size',[K M]}, mfilename, 'C0');
    end

    %% 3. Allocate histories and compute the fixed soft-thresholding scale
    % LambdaC penalizes change from the preceding coefficients. It changes
    % both the center and the threshold of each coefficient subproblem.
    max_iterations = options.MaxIterations;
    info.objective = zeros(max_iterations + 1, 1);
    info.nonzero_fraction = zeros(max_iterations + 1, 1);
    info.step_norm = zeros(max_iterations, 1);
    % Orthogonality gives ||X-D*C||_F = ||D'*X-C||_F. Reusing projected
    % avoids an extra D*C matrix multiplication for every objective value.
    info.objective(1) = 0.5*norm(projected - C, 'fro')^2 + beta*sum(abs(C(:)));
    info.nonzero_fraction(1) = nnz(C)/numel(C);
    tolerance_mode = strcmp(options.StopRule, 'tolerance');
    denominator = 1 + options.LambdaC;
    threshold = beta/denominator;

    %% 4. Alternate one dictionary update and one coefficient update
    for iteration = 1:max_iterations
        % X stays fixed. The term LambdaD*D draws the new basis toward the
        % previous basis while X*C' encodes the current patch-fitting term.
        % Exact orthogonal Procrustes solution of the dictionary subproblem.
        [U, ~, V] = svd(X*C' + options.LambdaD*D, 'econ');
        Dnew = U*V'; % Orthogonal by construction: Dnew'*Dnew = I.

        % Exact soft-threshold solution of the coefficient subproblem.
        projected = Dnew'*X;
        % Blend the newly projected patches with the preceding coefficients.
        center = (projected + options.LambdaC*C)/denominator;
        Cnew = sign(center).*max(abs(center) - threshold, 0); % Small coefficients become zero.

        % Combine the dictionary and coefficient changes into one absolute norm.
        step = hypot(norm(Dnew - D, 'fro'), norm(Cnew - C, 'fro'));
        % Commit both blocks before recording the objective and testing the stop.
        D = Dnew;
        C = Cnew;
        info.step_norm(iteration) = step;
        info.objective(iteration + 1) = ...
            0.5*norm(projected - C, 'fro')^2 + beta*sum(abs(C(:)));
        info.nonzero_fraction(iteration + 1) = nnz(C)/numel(C);
        % Fixed-count mode always completes its budget; tolerance mode may stop early.
        if tolerance_mode && step <= options.Tolerance
            break;
        end
    end

    %% 5. Report whether the requested inner accuracy was actually reached
    % Fixed-count mode may satisfy Tolerance, but its stopping reason remains
    % fixed_iterations because no tolerance test controlled that run.
    info.iterations = iteration;
    info.tolerance_met = step <= options.Tolerance;
    if ~tolerance_mode
        info.stop_reason = 'fixed_iterations';
    elseif info.tolerance_met
        info.stop_reason = 'tolerance';
    else
        info.stop_reason = 'max_iterations';
    end
    % Remove unused history rows if tolerance mode stopped before the cap.
    info.objective = info.objective(1:iteration + 1);
    info.nonzero_fraction = info.nonzero_fraction(1:iteration + 1);
    info.step_norm = info.step_norm(1:iteration);
    info.elapsed_seconds = toc(timer);
end

%% Settings and input checks (not part of alternating minimization)
function options = read_options(X,beta,options)
    validateattributes(X, {'double'}, {'real','finite','2d','nonempty','nonsparse'}, ...
        mfilename, 'X');
    validateattributes(beta, {'double'}, {'real','finite','scalar','nonnegative'}, ...
        mfilename, 'beta');
    validateattributes(options, {'struct'}, {'scalar'}, mfilename, 'options');

    defaults = struct('StopRule','fixed', 'MaxIterations',20, ...
        'Tolerance',1e-3, 'LambdaD',1, 'LambdaC',1, 'D0',[], 'C0',[]);
    names = fieldnames(options);
    % User fields override defaults; unrecognized fields raise an error.
    for j = 1:numel(names)
        if ~isfield(defaults, names{j})
            error('dictionary_reconstruction:UnknownOption', ...
                'Unknown option: %s.', names{j});
        end
        defaults.(names{j}) = options.(names{j});
    end
    options = defaults;
    options.StopRule = validatestring(options.StopRule, {'fixed','tolerance'}, ...
        mfilename, 'StopRule');
    validateattributes(options.MaxIterations, {'numeric'}, ...
        {'real','finite','scalar','integer','positive'}, mfilename, 'MaxIterations');
    validateattributes(options.Tolerance, {'numeric'}, ...
        {'real','finite','scalar','nonnegative'}, mfilename, 'Tolerance');
    validateattributes(options.LambdaD, {'numeric'}, ...
        {'real','finite','scalar','positive'}, mfilename, 'LambdaD');
    validateattributes(options.LambdaC, {'numeric'}, ...
        {'real','finite','scalar','positive'}, mfilename, 'LambdaC');
    % Avoid integer arithmetic if a numeric option was supplied as an integer.
    options.MaxIterations = double(options.MaxIterations);
    options.Tolerance = double(options.Tolerance);
    options.LambdaD = double(options.LambdaD);
    options.LambdaC = double(options.LambdaC);

end
