function [u,info] = blip_reconstruction(data,RF,TR,r,zeta,u0,options)
%BLIP_RECONSTRUCTION BLIP with an initial map, histories and bounded density.
%   Same signal-gradient/projection and step-halving rule as helper BLIP.m.
%   T1/T2 are selected from a Bloch dictionary; density is real and >=0.
%   Projection minimizes signal distance, including the density bounds.
%   Unlike the legacy routine, a rejected trial never replaces the output.
%   All pixels remain free: there is no support mask or spatial penalty.
%   To start from zero signal, use positive reference times and zero density.
%   MaxIterations counts ACCEPTED steps; MaxTrials includes rejected trials.
%   Times include dictionary construction, matching, and optional monitoring.
%   MATLAB R2021b; no additional toolbox is needed.
%
%   Algorithm reference:
%   M. Davies, G. Puy, P. Vandergheynst and Y. Wiaux,
%   "A Compressed Sensing Framework for Magnetic Resonance Fingerprinting,"
%   SIAM Journal on Imaging Sciences, 7(4), 2623-2656, 2014.
%   https://doi.org/10.1137/130947246
%   Code provenance: project adaptation of helper_functions/BLIP.m.
%   The legacy code's author, copyright holder and license are unverified.
%
%   Inputs: data is N-by-M-by-L complex k-space; u0 is N-by-M-by-3
%   physical [T1,T2,density]. RF, TR and zeta have L entries. RF gives flip
%   angles/phases, TR is in ms, and zeta selects acquired rows per frame.
%   Outputs: u has the same map ordering as u0. info records the initial
%   state and accepted updates, trial counts, and the reason for stopping.

    if nargin<7, options=struct(); end
    root=fileparts(mfilename('fullpath'));
    addpath(fullfile(root,'bloch_simulations'));

    %% 1. Read settings and check the experiment
    % Merge user options with defaults before allocating dictionary/history arrays.
    settings=read_options(options,data,RF,TR,r,zeta,u0);
    timer=tic;
    [n,m,L]=size(data); P=n*m;

    %% 2. Define the measured Fourier samples and optional noise threshold
    RF=RF(:).'; TR=TR(:).'; zeta=zeta(:).'; % Normalize vector shapes without conjugating RF.
    sampled=false(n,m,L);
    % Repeat the row selection over every Fourier column in each frame.
    for t=1:L, sampled(zeta(t):r:n,:,t)=true; end
    data(~sampled)=0; % Ignore unmeasured entries, also in the residual norm.
    threshold=NaN; % Comparisons with NaN are false, so [] disables noise stopping.
    if ~isempty(settings.NoiseLevel)
        % NoiseLevel is a norm over measured samples, not a standard deviation.
        assert(settings.NoiseLevel>=0 && isfinite(settings.NoiseLevel) && settings.DiscrepancyFactor>1, ...
            'Use a finite nonnegative noise norm and discrepancy factor >1.');
        threshold=settings.DiscrepancyFactor*settings.NoiseLevel;
    end

    %% 3. Simulate the discrete Bloch dictionary used for map projection
    % Enumerate every (T1,T2) pair; T2 changes fastest in this ordering.
    [t2,t1]=ndgrid(settings.T2Grid,settings.T1Grid);
    B=Blochmap_vectorized(RF,TR,t1(:),t2(:),0);
    % One row per (T1,T2) candidate; one column per acquisition frame.
    dict=reshape(B,numel(t1),L);
    norm_squared=sum(abs(dict).^2,2); % Needed for the optimal density of each atom.
    assert(all(norm_squared>0),'The RF sequence gives a zero dictionary signal.');
    dict_real=real(dict); dict_imag=imag(dict); % Reused by batched matching.

    %% 4. Initialize the signal, residual and diagnostic histories
    u=u0;
    signal=u(:,:,3).*Blochmap_vectorized(RF,TR,u(:,:,1),u(:,:,2),0);
    residual=data-forward(signal); % Measured minus predicted signal, zero off the sampling set.
    accepted=0; trials=0; mu=settings.InitialStep;
    % State histories include u0 at index 1; step histories have no initial entry.
    for name={'data_residual_norm','objective','time'}
        info.(name{1})=zeros(settings.MaxIterations+1,1);
    end
    info.monitor=zeros(settings.MaxIterations+1,0);
    info.step_size=zeros(settings.MaxIterations,1);
    info.backtracking=zeros(settings.MaxIterations,1);
    info.stop_reason='max_iterations';
    record_state();

    %% 5. Alternate a signal-gradient step and projection onto Bloch signals
    % There is no spatial prior: pixels interact through the Fourier data.
    % Commit a candidate only after the step-size inequality passes.
    for iteration=1:settings.MaxIterations
        if info.data_residual_norm(accepted+1)<=threshold
            info.stop_reason='discrepancy'; break;
        end
        gradient_step=sqrt(P)*ifft2(residual); % Adjoint of unitary sampled FFT.
        success=false;
        for trial=1:settings.MaxBacktracking
            if trials>=settings.MaxTrials, break; end
            trials=trials+1;
            % Match each pixel to a dictionary atom and its best bounded density.
            [candidate,candidate_signal]=project(signal+mu*gradient_step);
            change=candidate_signal-signal;
            % Compare image/time-space change with its change on acquired samples.
            signal_change=sum(abs(change(:)).^2);
            measured_change=forward(change);
            data_change=sum(abs(measured_change(:)).^2);
            % Cross-multiplied legacy mu <= kappa*||d||^2/||A*d||^2.
            % This form also handles a zero step without dividing by zero.
            if mu*data_change<=settings.Kappa*signal_change+eps*max(1,signal_change)
                success=true; break;
            end
            mu=mu/settings.BacktrackingFactor; % Retry with a smaller signal step.
        end
        if ~success
            info.stop_reason='backtracking_failed';
            if trials>=settings.MaxTrials, info.stop_reason='max_trials'; end
            break; % Retain the preceding accepted maps and signal.
        end
        % Replace the current maps only after the step-size check succeeds.
        relative_step=sqrt(signal_change)/(1+norm(signal(:)));
        u=candidate; signal=candidate_signal;
        residual=residual-measured_change; % Reuse FFT from the acceptance check.
        accepted=accepted+1;
        info.step_size(accepted)=mu; info.backtracking(accepted)=trial-1;
        record_state();
        if settings.Verbose
            fprintf('%3d BLIP residual=%.6g step=%.2g trials=%d\n', ...
                accepted,info.data_residual_norm(accepted+1),mu,trial);
        end
        if info.data_residual_norm(accepted+1)<=threshold
            info.stop_reason='discrepancy'; break;
        elseif signal_change==0 || (settings.StepTolerance>0 && relative_step<=settings.StepTolerance)
            info.stop_reason='small_signal_step'; break;
        end
    end

    %% 6. Return only histories belonging to accepted updates
    for name={'data_residual_norm','objective','time','monitor'}
        info.(name{1})=info.(name{1})(1:accepted+1,:);
    end
    info.step_size=info.step_size(1:accepted); info.backtracking=info.backtracking(1:accepted);
    info.iterations=accepted; info.total_trials=trials;
    info.elapsed_seconds=toc(timer); info.options=settings;

    %% Helper: unitary FFT restricted to the measured rows
    function measured=forward(values)
        measured=fft2(values)/sqrt(P);
        measured(~sampled)=0;
    end

    %% Helper: closest bounded Bloch signal for each pixel
    function [maps,matched]=project(values)
        values=reshape(values,P,L).'; % One temporal signal per column.
        index=zeros(P,1); rho=zeros(P,1);
        for first=1:settings.BatchSize:P
            % Batching limits the dictionary-by-pixel correlation matrix in memory.
            pixels=first:min(P,first+settings.BatchSize-1);
            % Real inner products permit real, nonnegative density estimates.
            inner=dict_real*real(values(:,pixels))+dict_imag*imag(values(:,pixels));
            % For each atom, minimize over density analytically, then clip to bounds.
            amplitudes=min(settings.DensityBounds(2),max(settings.DensityBounds(1),inner./norm_squared));
            % ||v-rho*b||^2 = ||v||^2 - (2*rho*<b,v>-rho^2*||b||^2).
            score=2*amplitudes.*inner-amplitudes.^2.*norm_squared;
            [~,chosen]=max(score,[],1); % Maximum score gives minimum signal distance.
            entries=sub2ind(size(amplitudes),chosen,1:numel(pixels));
            % Keep each winning atom together with that atom's optimal density.
            index(pixels)=chosen; rho(pixels)=amplitudes(entries);
        end
        maps=cat(3,reshape(t1(index),n,m),reshape(t2(index),n,m),reshape(rho,n,m));
        matched=reshape(rho.*dict(index,:),n,m,L);
        % T1/T2 have no signal meaning where rho=0; keep their previous values.
        for channel=1:2
            old=u(:,:,channel); current=maps(:,:,channel);
            current(rho==0)=old(rho==0); maps(:,:,channel)=current;
        end
    end

    %% Helper: record measured-data error and optional monitoring
    function record_state()
        % These are diagnostics of the current accepted state, not update formulas.
        k=accepted+1;
        info.data_residual_norm(k)=norm(residual(:));
        info.objective(k)=.5*info.data_residual_norm(k)^2;
        if ~isempty(settings.MonitorFcn)
            metrics=settings.MonitorFcn(u); info.monitor(k,1:numel(metrics))=metrics(:).';
        end
        info.time(k)=toc(timer);
    end
end

%% Settings and input checks (not part of the BLIP iteration)
function settings=read_options(options,data,RF,TR,r,zeta,u0)
    % assert(condition,message) stops with that message when condition is false.
    defaults=struct('T1Grid',[100:100:2000 2250:250:5500], ...
        'T2Grid',[10:5:100 125:25:600],'DensityBounds',[0 110], ...
        'MaxIterations',40,'MaxTrials',Inf,'MaxBacktracking',30, ...
        'InitialStep',r,'BacktrackingFactor',2,'Kappa',.99, ...
        'StepTolerance',0,'NoiseLevel',[],'DiscrepancyFactor',1.1, ...
        'BatchSize',512,'MonitorFcn',[],'Verbose',true);
    names=fieldnames(options);
    % Reject unknown options so misspellings cannot silently use a default.
    for j=1:numel(names)
        assert(isfield(defaults,names{j}),'Unknown BLIP option: %s.',names{j});
        defaults.(names{j})=options.(names{j});
    end
    settings=defaults;
    [n,m,L]=size(data);
    % Check dimensions, sampling, physical feasibility and iteration budgets.
    assert(isequal(size(u0),[n m 3]) && numel(RF)==L && numel(TR)==L, ...
        'Use N-by-M-by-3 initial maps and one RF/TR value per frame.');
    assert(r>=1 && r==fix(r) && mod(n,r)==0 && numel(zeta)==L && ...
        all(zeta>=1 & zeta<=r & zeta==fix(zeta)), 'Invalid Fourier sampling.');
    initial_times=u0(:,:,1:2);
    assert(all(isfinite(u0(:))) && all(initial_times(:)>0) && ...
        settings.DensityBounds(1)>=0 && settings.DensityBounds(2)>=settings.DensityBounds(1), ...
        'Use positive times and nonnegative density bounds.');
    density=u0(:,:,3);
    assert(all(density(:)>=settings.DensityBounds(1) & density(:)<=settings.DensityBounds(2)), ...
        'Initial density must satisfy bounds.');
    assert(settings.InitialStep>0 && settings.BacktrackingFactor>1 && settings.Kappa>0 && settings.Kappa<1, ...
        'Use positive step, backtracking factor >1 and 0<Kappa<1.');
    assert(settings.BatchSize>=1 && settings.BatchSize==fix(settings.BatchSize) && ...
        settings.MaxIterations>=0 && settings.MaxIterations==fix(settings.MaxIterations) && ...
        settings.MaxBacktracking>=1 && settings.MaxBacktracking==fix(settings.MaxBacktracking) && ...
        settings.MaxTrials>=0 && settings.MaxTrials==fix(settings.MaxTrials), 'Invalid iteration or batch limit.');
    assert(~isempty(settings.T1Grid) && ~isempty(settings.T2Grid) && ...
        all(isfinite(settings.T1Grid) & settings.T1Grid>0) && all(isfinite(settings.T2Grid) & settings.T2Grid>0), ...
        'Dictionary grids must contain positive finite times.');
end
