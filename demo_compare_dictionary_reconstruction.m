%% Compare nested and exactly one-step dictionary learning
% Run in MATLAB R2021b + Optimization Toolbox. Optional artificial maps are saved below.
% Both methods use the same measurements, initial maps, D0/C0, objective,
% map solver and outer stopping tolerances. Only the inner update rule differs.
% Optional saved 128/256 presets: see the selector after section 1 below.
% Workflow: synthesize data -> create one common u0/D0/C0 -> run both rules
% through nested_dictionary_reconstruction -> compare final maps and histories.
% Clear previous variables/figures, then seed random noise for reproducible runs.
clearvars;
close all;
rng(1);
% Locate data and helpers relative to this script's own folder.
root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'helper_functions'),fullfile(root,'bloch_simulations'));

%% Input: change the file and physical ranges here
% MAT file with equally sized 2-D arrays T1, T2 (ms) and nonnegative Density.
% Density=0 marks background. Save your maps with save('my_maps.mat','T1','T2','Density').
% These are reference maps for simulated measurements; optional rescaling is selected below.
input_file=fullfile(root,'data','PWbrain256.mat'); % Or the full path to your own MAT file.
% OPTIONAL artificial input: uncomment ALL THREE lines below to use it.
% Higher complexity adds more/smaller regions and finer texture. The MAT file is regenerated each run.
% [T1,T2,Density]=create_artificial_ground_truth(10,256); % Complexity, image size.
% input_file=fullfile(root,'data','artificial_maps.mat');
% save(input_file,'T1','T2','Density');
% DOWNSAMPLING a 256-by-256 input: 1 -> 256-by-256, 2 -> 128-by-128, 4 -> 64-by-64.
pixel_stride=2;              % Keep every pixel_stride-th row/column; no interpolation.
% Optional direct division is shown after loading the maps below.
lower_bounds=[100 10 0];     % Physical [T1(ms), T2(ms), Density]; positive times, zero density.
upper_bounds=[5500 600 110]; % Edit for the range expected in your own maps.
scales=[5000 500 100]; reference_times=[1000 100];
% The solvers learn patches of (u-[reference_times,0])./scales. These choices
% therefore set both numerical conditioning and the units of regularization.
% Keep reference_times inside the time bounds; review prior weights if scales change.
initial_values=[2000 200 50]; % Uniform foreground start, inside all three bounds.
% If you divide a map below, adjust its bounds, scales, reference_times
% and initial_values manually; these settings are not changed automatically.

%% 1. Shared acquisition and solver settings
L=100; r=8; TR_ms=12;         % r must divide the reduced row count; r=1 works for any size.
noise_std=.8;                % Standard deviation of each real/imaginary noise component.
patch_size=[8 8];
alpha=.1;                   % Spatial penalty on shifted, scaled maps.
lambda=[4.5 45 45]*3.5;     % Dictionary weights [T1,T2,density].
beta=[.0045 .0045 .0045]*5; % L1 weights inside Lambda*(fit/2 + Beta*|C|_1).
max_iterations=100;         % Same outer budget; elapsed times can differ.
gauss_newton='approximate'; % 'exact' uses the sampled Fourier Hessian.
step_norm='H01';            % 'H01': gradients; 'H1': values plus gradients.
inner_max_iterations=2000;  % Cap for NESTED only; one-step always uses one pair.
inner_tolerance_scale=.007;
% Nested accuracy decreases as scale*||(D0,C0)||/(k+1)^inner_decay.
% The one-step rule performs exactly one D/C pair regardless of this tolerance.
inner_decay=.65;
qp_optimality_tolerance=1e-8;
% Separate the QP's stationarity test, objective-change test and inner PCG accuracy.
qp_function_tolerance=0;
qp_pcg_tolerance=1e-6;

%% Optional saved qMRI parameters
% [] keeps this demo's manual settings above. Uncomment ONE preset to load it.
% Presets assume a 256-by-256 input; they override acquisition, bounds/scales
% and solver settings. The input file and uniform initialization stay as above.
parameters=[];
% parameters=reconstruction_parameters_128(); % Saved 128-by-128 baseline.
% parameters=reconstruction_parameters_256(); % Independent 256-by-256 copy to tune.
% Put one-run changes to parameters HERE, before the assignments below.
if ~isempty(parameters)
    pixel_stride=parameters.pixel_stride;
    lower_bounds=parameters.lower_bounds; upper_bounds=parameters.upper_bounds;
    reference_times=parameters.reference_times; scales=parameters.scales;
    L=parameters.L; r=parameters.r; TR_ms=parameters.TR_ms;
    noise_std=parameters.noise_std;
    rng(parameters.seed);
    % This controlled two-method demo uses nested_options for BOTH methods;
    % only InnerStopRule differs. Use the four-method demo for separate options.
    patch_size=parameters.nested_options.PatchSize;
    beta=parameters.nested_options.Beta;
end

%% Plot design: same defaults as demo_compare_all_reconstructions
% These options affect display copies and PDF exports only.
design=reconstruction_plot_design(lower_bounds,upper_bounds,[4 3]);
design.ShowFigures=true; design.ShowRowFigures=false;
design.ZeroBackground=true; design.RowShowTitles=false;
design.ProfileRow=[]; % [] chooses the middle row.
design.ZoomBox=[]; % [first_row last_row first_col last_col]; [] selects a body region.
% MapLimits, ErrorLimits, MapColors, ErrorColors and sizes are editable in design.
output_folder=fullfile(root,'data','demo_reports',['comparison_dictionary_' datestr(now,'yyyymmdd_HHMMSS')]);

%% 2. Generate one common data set
% Load physical reference maps. Positive background T1/T2 values keep the
% Bloch model valid; zero background density still makes its signal zero.
maps=load(input_file,'T1','T2','Density');
% Convert numeric storage classes to double; all maps must describe the same grid.
for name={'T1','T2','Density'}
    validateattributes(maps.(name{1}),{'numeric'},{'real','finite','nonempty','2d'});
    maps.(name{1})=double(maps.(name{1}));
end
assert(isequal(size(maps.T1),size(maps.T2),size(maps.Density)), ...
    'T1, T2 and Density must have the same size.');
assert(all(maps.Density(:)>=0),'Density must be nonnegative.');
% OPTIONAL: uncomment this line to use T1/10 as the new ground truth.
% Adjust the settings above manually. Use the same form for T2 or Density if wanted.
% maps.T1 = maps.T1 / 10;
outside=maps.Density==0;
maps.T1(outside)=reference_times(1); maps.T2(outside)=reference_times(2);
% Sample the same spatial grid in each map, then stack channels as T1,T2,rho.
truth=cat(3,maps.T1(1:pixel_stride:end,1:pixel_stride:end), ...
    maps.T2(1:pixel_stride:end,1:pixel_stride:end), ...
    maps.Density(1:pixel_stride:end,1:pixel_stride:end));
mask=truth(:,:,3)>0;
[n,m,~]=size(truth);
% Reject a mismatched input range before simulating measurements.
assert(any(mask(:)),'The reduced input must contain positive Density values.');
for j=1:3
    values=truth(:,:,j);
    assert(all(values(mask)>=lower_bounds(j) & values(mask)<=upper_bounds(j)), ...
        'Input map %d is outside the physical bounds. Edit lower_bounds/upper_bounds.',j);
end
assert(mod(n,r)==0,'r must divide the reduced image row count; use r=1 for any size.');
% Load only the RF pulses. TR is set separately and is shared by simulation
% and reconstruction; zeta cycles the starting row of the sampling pattern.
sequence=load(fullfile(root,'data','sequence_parameter_middle.mat'),'RF');
RF=sequence.RF(1:2:2*L); TR=TR_ms*ones(1,L); zeta=mod(0:L-1,r)+1;
% Simulate a unit-density complex response for every pixel and time frame.
B=Blochmap_vectorized(RF,TR,truth(:,:,1),truth(:,:,2),0);
% Allocate full k-space storage; only measured rows receive signal and noise.
data=complex(zeros(n,m,L));
% Multiply by density before the unitary FFT. Simulate each time frame separately.
for t=1:L
    signal=fft2(truth(:,:,3).*B(:,:,t))/sqrt(n*m);
    % Keep every r-th Fourier row. Independent real/imaginary Gaussian draws
    % give expected squared complex noise magnitude 2*noise_std^2 per sample.
    rows=zeta(t):r:n;
    data(rows,:,t)=signal(rows,:)+noise_std*(randn(numel(rows),m)+1i*randn(numel(rows),m));
end

%% 3. Common initialization from u0, without learning from the true maps
offset=reshape([reference_times 0],1,1,3);
% Start the background at the reference state and the foreground at uniform
% initial_values. No foreground ground-truth parameter values are copied.
u0=repmat(offset,n,m,1);
for j=1:3
    initial=u0(:,:,j); initial(mask)=initial_values(j); u0(:,:,j)=initial;
end
x0=(u0-offset)./reshape(scales,1,1,3);
% Build shared dictionary initializations from these starting deviations.
[R,number_of_patches]=create_R_matrix([n m],patch_size,1);
patch_area=prod(patch_size);
D0=zeros(patch_area,patch_area,3); C0=zeros(patch_area,number_of_patches,3);
for j=1:3
    X=create_data_matrix_fast(x0(:,:,j),patch_size,R);
    [D0(:,:,j),~,~]=svd(X*X','econ');
    % Project initial patches onto this orthogonal PCA basis and soft-threshold.
    C0(:,:,j)=prox_c(beta(j),D0(:,:,j)'*X,'soft');
end

%% 4. Run the two methods with the same settings except InnerStopRule
% MonitorFcn evaluates errors only; truth never enters the optimization.
options=struct('Scales',scales,'ReferenceTimes',reference_times,'H',1, ...
    'Alpha',alpha,'Lambda',lambda,'Beta',beta,'PatchSize',patch_size, ...
    'StepNorm',step_norm,'GaussNewton',gauss_newton,'D0',D0,'C0',C0, ...
    'LowerBounds',lower_bounds,'UpperBounds',upper_bounds,'BodyMask',mask, ...
    'MaxIterations',max_iterations,'MapTolerance',1e-6,'DictionaryTolerance',1e-3, ...
    'InnerMaxIterations',inner_max_iterations,'InnerToleranceScale',inner_tolerance_scale, ...
    'InnerDecay',inner_decay,'QPOptimalityTolerance',qp_optimality_tolerance, ...
    'QPFunctionTolerance',qp_function_tolerance,'QPPCGTolerance',qp_pcg_tolerance, ...
    'QPRequireOptimality',true,'QPMaxIterations',300,'Verbose',true, ...
    'MonitorFcn',@(u) relative_body_error(u,truth,mask));
if ~isempty(parameters)
    % Copy all saved solver controls; retain the runtime maps/mask above.
    for field=fieldnames(parameters.nested_options).'
        name=field{1}; options.(name)=parameters.nested_options.(name);
    end
    options.H=parameters.h;
end
rules={'tolerance','one-step'}; method_names={'Nested DL','One-step DL'};
% Both calls receive identical initial arrays. MATLAB passes them by value,
% so learning in the first call cannot alter the second call's starting D/C.
reconstructions=cell(1,2); dictionaries=cell(1,2); coefficients=cell(1,2); histories=cell(1,2);
for method=1:2
    options.InnerStopRule=rules{method};
    % Switch only how far dictionary learning runs inside each outer map update.
    fprintf('\n%s: %dx%d, L=%d, r=%d, %s GN, %s step norm\n', ...
        method_names{method},n,m,L,r,options.GaussNewton,options.StepNorm);
    [reconstructions{method},dictionaries{method},coefficients{method},histories{method}]= ...
        nested_dictionary_reconstruction(data,RF,TR,r,zeta,u0,options);
    state=histories{method};
    fprintf('Stopped by %s after %d outer steps, %.2f seconds.\n', ...
        state.stop_reason,state.iterations,state.elapsed_seconds);
    if strcmp(state.stop_reason,'inner_tolerance_not_met')
        fprintf('Last inner changes: [%g %g %g]; requested tolerance: %g.\n', ...
            state.last_inner_step_norm,state.last_inner_tolerance);
    end
end
% Each timer includes solver setup and diagnostic work, but excludes common
% data generation/D0/C0 construction and plotting. These are single-run
% timings; use repeated runs with alternating order for formal benchmarks.

%% 5. Final errors, objective, runtime and the actual stopping reason
outer_steps=zeros(2,1); seconds=zeros(2,1); final_objective=zeros(2,1);
relative_error=zeros(2,3); stop_reason=cell(2,1);
for method=1:2
    state=histories{method};
    % State histories contain the start and every completed outer update;
    % end selects the actual final state even after early termination.
    outer_steps(method)=state.iterations; seconds(method)=state.elapsed_seconds;
    final_objective(method)=state.objective(end);
    relative_error(method,:)=state.monitor(end,:); stop_reason{method}=state.stop_reason;
end
comparison=table(outer_steps,seconds,final_objective,relative_error(:,1), ...
    relative_error(:,2),relative_error(:,3),stop_reason, ...
    'VariableNames',{'OuterSteps','Seconds','Objective','T1Error_percent', ...
    'T2Error_percent','DensityError_percent','StopReason'},'RowNames',method_names);
disp(comparison);
fprintf('Initial body errors [T1 T2 density]: [%.2f %.2f %.2f] percent.\n', ...
    relative_body_error(u0,truth,mask));

%% 6. Shared map/error reports, zoom insets, compact rows and profiles
export_reconstruction_report(truth,reconstructions,u0,mask,method_names,design,output_folder);
labels=design.MapNames;

%% 7. Compare optimization progress by outer steps AND elapsed time
% Nested steps can do more dictionary work, so iteration count and elapsed
% time answer different questions. max(...,eps) only keeps log plots finite.
colors=design.MethodColors;
fig=reconstruction_report_figure('Dictionary comparison: convergence',design);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile; hold on;
for method=1:2
    state=histories{method};
    plot(0:state.iterations,max(state.objective,eps),'Color',colors(method,:),'LineStyle',design.LineStyles{method},'LineWidth',design.LineWidth);
end
set(gca,'YScale','log'); grid on; xlabel('Outer step'); title('Full objective');
legend(method_names,'Location','best');
nexttile; hold on;
for method=1:2
    state=histories{method};
    plot(state.time,max(state.objective,eps),'Color',colors(method,:),'LineStyle',design.LineStyles{method},'LineWidth',design.LineWidth);
end
set(gca,'YScale','log'); grid on; xlabel('Elapsed seconds'); title('Full objective versus time');
nexttile; hold on;
for method=1:2
    state=histories{method};
    plot(0:state.iterations,max(state.projected_gradient,eps),'Color',colors(method,:),'LineStyle',design.LineStyles{method},'LineWidth',design.LineWidth);
end
set(gca,'YScale','log'); grid on; xlabel('Outer step'); title('Map projected gradient (D,C fixed)');
nexttile; hold on;
for method=1:2
    state=histories{method};
    stairs(1:state.iterations,sum(state.inner_iterations,2),'Color',colors(method,:),'LineStyle',design.LineStyles{method},'LineWidth',design.LineWidth);
end
grid on; xlabel('Completed outer step'); ylabel('D/C update pairs, summed over maps');
title('Dictionary work per outer step'); % One-step uses 3 pairs: one for each active map.

export_reconstruction_pdf(fig,fullfile(output_folder,'convergence.pdf'),design);

%% 8. Reconstruction errors versus iterations and time
% The MonitorFcn saved three scalar body errors per state; there is no need
% to store a full image after every iteration to draw these curves.
fig=reconstruction_report_figure('Dictionary comparison: errors',design);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
for row=1:2
    for j=1:3
        nexttile; hold on;
        for method=1:2
            state=histories{method};
            horizontal=0:state.iterations;
            if row==2, horizontal=state.time; end
            plot(horizontal,state.monitor(:,j),'Color',colors(method,:),'LineStyle',design.LineStyles{method},'LineWidth',design.LineWidth);
        end
        grid on; ylabel('Relative error inside body (%)'); title(labels{j});
        if row==1, xlabel('Outer step'); else, xlabel('Elapsed seconds'); end
    end
end
legend(method_names,'Location','best');

export_reconstruction_pdf(fig,fullfile(output_folder,'convergence.pdf'),design,true);

function errors=relative_body_error(u,truth,mask)
% Three scalar errors; storing full parameter maps at every step is unnecessary.
% Logical mask indexing restricts both the error norm and reference norm to tissue.
    errors=zeros(1,3);
    for j=1:3
        estimate=u(:,:,j); target=truth(:,:,j);
        errors(j)=100*norm(estimate(mask)-target(mask))/norm(target(mask));
    end
end
