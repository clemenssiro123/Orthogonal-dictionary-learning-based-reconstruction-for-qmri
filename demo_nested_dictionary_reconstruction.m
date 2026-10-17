%% Nested dictionary-learning qMRI reconstruction
% MATLAB R2021b + Optimization Toolbox. Run this script directly.
% The known body support is used for initialization, errors and display only.
% The solver learns its dictionaries from the current reconstruction.
% Optional saved 128/256 presets: see the selector after section 1 below.
% Workflow: choose acquisition/prior settings -> simulate measurements ->
% construct a uniform start -> alternate D/C learning and bounded map updates
% -> inspect map errors, objective terms, inner accuracy and spatial profiles.
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
scales=[5000 500 100];
% One scale per physical channel. The reference times are the origin of
% scaled deviations; zero deviation does not mean zero physical T1 or T2.
% Changing scales changes the prior, so review Alpha/Lambda/Beta if units change.
reference_times=[1000 100];   % Background times, inside the time bounds above.
initial_values=[2000 200 50]; % Uniform foreground start, inside all three bounds.
% If you divide a map below, adjust its bounds, scales, reference_times
% and initial_values manually; these settings are not changed automatically.

%% 1. Acquisition and solver settings
L=100; r=8; TR_ms=12;         % r must divide the reduced row count; r=1 works for any size.
noise_std=.8;                % Standard deviation of each complex-noise component.
patch_size=[8 8];            % Small patches for this reduced-size demonstration.
alpha=0.1;                    % Spatial gradient penalty on the scaled deviations.
lambda=[4.5 45 45]*3.5;          % Patch-fit weights, ordered [T1,T2,density].
beta=[.0045 .0045 .0045]*5;    % L1 weights INSIDE lambda*(fit/2 + beta*|C|_1).
max_iterations=100;
% This limits complete outer updates, not total D/C pairs or QP attempts.
gauss_newton='approximate';        % 'approximate' uses the paper's 1/r local Gram.
step_norm='H01';              % 'H01': gradients only; 'H1': values plus gradients.
inner_stop_rule='tolerance'; % 'fixed' with inner_max_iterations=1 gives one-step learning.
inner_max_iterations=500;
% A tolerance-controlled inner solve that reaches this cap without meeting
% eta stops the outer run and preserves the preceding complete iterate.
inner_tolerance_scale=.007;   % Smaller than 1 to demonstrate multiple inner updates.
inner_decay=.65;
% eta_k = inner_tolerance_scale*||(D0,C0)||/(k+1)^inner_decay.
qp_optimality_tolerance=1e-8;
qp_function_tolerance=0;
qp_pcg_tolerance=1e-6;
% These weights are starting values, not a reproduction of the manuscript's
% experiments: we use physical times, shifted coordinates and a smaller grid.

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
end

%% Plot design: same defaults as demo_compare_all_reconstructions
% These options affect display copies and PDF exports only.
design=reconstruction_plot_design(lower_bounds,upper_bounds,4);
design.ShowFigures=true; design.ShowRowFigures=false;
design.ZeroBackground=true; design.RowShowTitles=false;
design.ProfileRow=[]; % [] chooses the middle row.
design.ZoomBox=[]; % [first_row last_row first_col last_col]; [] selects a body region.
% MapLimits, ErrorLimits, MapColors, ErrorColors and sizes are editable in design.
output_folder=fullfile(root,'data','demo_reports',['nested_dictionary_' datestr(now,'yyyymmdd_HHMMSS')]);

%% 2. Physical ground truth and sampled noisy data
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

%% 3. Initialize physical maps and reconstruct
u0=repmat(reshape([reference_times 0],1,1,3),n,m,1);
% Use reference background values and fill tissue with uniform initial_values.
for j=1:3
    initial=u0(:,:,j); initial(mask)=initial_values(j); u0(:,:,j)=initial;
end
% Pack acquisition coordinates, regularization and QP/outer stopping controls.
options=struct('Scales',scales,'ReferenceTimes',reference_times, ...
    'Alpha',alpha,'Lambda',lambda,'Beta',beta,'PatchSize',patch_size, ...
    'H',1,'StepNorm',step_norm,'GaussNewton',gauss_newton, ...
    'LowerBounds',lower_bounds,'UpperBounds',upper_bounds,'BodyMask',mask, ...
    'MaxIterations',max_iterations,'MapTolerance',1e-6,'DictionaryTolerance',1e-3, ...
    'InnerStopRule',inner_stop_rule,'InnerMaxIterations',inner_max_iterations, ...
    'InnerToleranceScale',inner_tolerance_scale,'InnerDecay',inner_decay, ...
    'QPOptimalityTolerance',qp_optimality_tolerance, ...
    'QPFunctionTolerance',qp_function_tolerance,'QPPCGTolerance',qp_pcg_tolerance, ...
    'QPRequireOptimality',true,'QPMaxIterations',300,'Verbose',true);
if ~isempty(parameters)
    % Copy all saved solver controls; retain the runtime maps/mask above.
    for field=fieldnames(parameters.nested_options).'
        name=field{1}; options.(name)=parameters.nested_options.(name);
    end
    options.H=parameters.h;
end
[u,D,C,info]=nested_dictionary_reconstruction(data,RF,TR,r,zeta,u0,options);
% u contains physical maps; D/C describe patches of the solver's scaled maps.
% info records only complete outer updates, plus diagnostics of a failed trial.
fprintf('Stopped by %s after %d complete outer steps, %.2f seconds.\n', ...
    info.stop_reason,info.iterations,info.elapsed_seconds);

%% 4. Errors use the unmasked reconstruction inside the body
% Compare initial and final estimates using 100*||error(mask)||/||truth(mask)||.
relative_error=zeros(2,3);
for j=1:3
    target=truth(:,:,j); initial=u0(:,:,j); result=u(:,:,j);
    relative_error(:,j)=100*[norm(initial(mask)-target(mask)); ...
        norm(result(mask)-target(mask))]/norm(target(mask));
end
disp(array2table(relative_error,'VariableNames',{'T1','T2','Density'}, ...
    'RowNames',{'Initial','Nested_dictionary'}));

%% 5. Shared map/error reports, zoom insets, compact rows and profiles
method_names={'Nested DL'};
export_reconstruction_report(truth,{u},u0,mask,method_names,design,output_folder);

%% 6. Outer convergence and accuracy of the inner dictionary solves
% State histories include iteration zero; step/inner histories start at one.
% The final inner D/C change is different from the total outer D/C change.
% The inner-tolerance plot compares each channel's LAST inner change with eta_k.
k=0:info.iterations;
fig=reconstruction_report_figure('Nested dictionary convergence',design);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
nexttile;
colororder(gca,lines(5));
semilogy(k,max([info.objective info.data info.patch info.spatial info.sparsity],eps),'LineWidth',design.LineWidth);
grid on; xlabel('Outer step'); title('Objective contributions');
legend('Total','Data','Patch fit','Spatial','L1','Location','best');
nexttile;
semilogy(k,max(info.projected_gradient,eps),'LineWidth',design.LineWidth);
grid on; xlabel('Outer step'); title('Map projected gradient (D,C fixed)');
nexttile;
colororder(gca,lines(2));
semilogy(1:info.iterations,max([info.step_norm info.dictionary_step_norm],eps),'LineWidth',design.LineWidth);
grid on; xlabel('Outer step'); title('Absolute outer changes'); legend('Map','D,C');
nexttile;
colororder(gca,lines(3));
plot(1:info.iterations,info.inner_iterations,'LineWidth',design.LineWidth);
grid on; xlabel('Outer step'); title('Inner D/C update pairs'); legend('T1','T2','Density');
nexttile;
colororder(gca,lines(3));
semilogy(1:info.iterations,max(info.inner_step_norm,eps),'LineWidth',design.LineWidth); hold on;
semilogy(1:info.iterations,info.inner_tolerance,'k--','LineWidth',design.LineWidth);
grid on; xlabel('Outer step'); title('Last inner change versus eta_k');
legend('T1','T2','Density','eta_k','Location','best');
nexttile;
plot(k,info.background_rms(:,3),'LineWidth',design.LineWidth);
grid on; xlabel('Outer step'); title('Background RMS density');

export_reconstruction_pdf(fig,fullfile(output_folder,'convergence.pdf'),design);
