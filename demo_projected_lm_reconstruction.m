%% Projected Levenberg-Marquardt qMRI reconstruction (MATLAB R2021b)
% Uses the Bloch signal model followed by sampled unitary Fourier transforms.
% DATA-ONLY objective. L2 damping penalizes PHYSICAL increments, as in the
% pasted algorithm. Each update is an unconstrained solve followed by clipping.
% Damping follows a schedule; objective increases are retained.
% Known body support is used for initialization, diagnostics and display.
% Optional saved 128/256 presets: see the selector after section 1 below.
% Workflow: generate one noisy data set -> initialize physical maps -> solve
% damped normal equations and clip to bounds -> inspect fit, errors and steps.
% Clear previous variables/figures, then seed random noise for reproducible runs.
clearvars;
close all;
rng(1);
% Locate data and helpers relative to this script's own folder.
root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'plot_routines'));
addpath(fullfile(root,'bloch_simulations'),fullfile(root,'helper_functions'));

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
pixel_stride=1;           % Keep every pixel_stride-th row/column; no interpolation.
% Optional direct division is shown after loading the maps below.
lower_bounds=[100 10 0];  % Physical [T1(ms), T2(ms), Density]; positive times, zero density.
upper_bounds=[5500 600 110]; % Edit for the range expected in your own maps.
scales=[5000 500 100];
% Scale each channel to dimensionless variables. The solver's damping metric
% compensates for these scales so damping is defined on physical increments.
reference_times=[1000 100]; % Background times, inside the time bounds above.
initial_values=[2000 200 50]; % Uniform foreground start, inside all three bounds.
% If you divide a map below, adjust its bounds, scales, reference_times
% and initial_values manually; these settings are not changed automatically.

%% 1. Acquisition and solver settings
L=100;
r=32;                      % Must divide the reduced row count; r=1 works for any size.
noise_std=1.0;             % Standard deviation of EACH real/imaginary part.
TR_ms=12;
max_iterations=14;
gauss_newton='approximate'; % Sparse direct solve; 'exact' uses sampled FFT + PCG.
initial_damping=1;        % Damping before the first update (k=0).
damping_decay=.6;         % lambda_k=max(floor,lambda_0*decay^k,weight*residual).
residual_damping_weight=1e-12;
min_damping=1e-8;
linear_tolerance=1e-8;    % Checked relative residual BEFORE projection.
linear_max_iterations=500; % PCG limit for exact GN.
stationarity_tolerance=1e-6;
% This controls x-clip(x-gradient), while linear_tolerance controls the
% accuracy of each linear system. They are different stopping quantities.
step_tolerance=0;         % Disable stopping solely because a step is small.
use_discrepancy=false;     % Stop at the first measured residual <= tau*delta.
discrepancy_factor=1.002;   % tau>1, independently adjustable.
noise_level=[];          % delta: noise NORM, not sigma; [] estimates it below.
% With no spatial prior there are no exterior Dirichlet conditions. The
% reference times are just the coordinate origin and background initialization.
% See Dong et al. (2019), sec. 4.4: https://arxiv.org/pdf/1903.11979

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
    use_discrepancy=~isempty(parameters.lm_options.NoiseLevel);
    noise_level=parameters.lm_options.NoiseLevel;
end

%% Plot design: same defaults as demo_compare_all_reconstructions
% These options affect display copies and PDF exports only.
design=reconstruction_plot_design(lower_bounds,upper_bounds,2);
design.ShowFigures=true; design.ShowRowFigures=false;
design.ZeroBackground=true; design.RowShowTitles=false;
design.ProfileRow=[]; % [] chooses the middle row.
design.ZoomBox=[]; % [first_row last_row first_col last_col]; [] selects a body region.
% MapLimits, ErrorLimits, MapColors, ErrorColors and sizes are editable in design.
output_folder=fullfile(root,'data','demo_reports',['projected_lm_' datestr(now,'yyyymmdd_HHMMSS')]);

%% 2. Physical truth and sampled Fourier data
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
source_mask=maps.Density>0;
maps.T1(~source_mask)=reference_times(1);
maps.T2(~source_mask)=reference_times(2);
maps.Density(~source_mask)=0;
% Sample the same spatial grid in each map, then stack channels as T1,T2,rho.
truth=cat(3,maps.T1(1:pixel_stride:end,1:pixel_stride:end), ...
    maps.T2(1:pixel_stride:end,1:pixel_stride:end), ...
    maps.Density(1:pixel_stride:end,1:pixel_stride:end));
mask=truth(:,:,3)>0;
background_values=[reference_times 0];
[n,m,~]=size(truth);
% Reject a mismatched input range before simulating measurements.
assert(any(mask(:)),'The reduced input must contain positive Density values.');
for j=1:3
    values=truth(:,:,j);
    assert(all(values(mask)>=lower_bounds(j) & values(mask)<=upper_bounds(j)), ...
        'Input map %d is outside the physical bounds. Edit lower_bounds/upper_bounds.',j);
end
assert(mod(n,r)==0,'The undersampling rate r must divide the row count.');
% Load only the RF pulses. TR is set separately and is shared by simulation
% and reconstruction; zeta cycles the starting row of the sampling pattern.
sequence=load(fullfile(root,'data','sequence_parameter_middle.mat'),'RF');
assert(2*L-1<=numel(sequence.RF),'L exceeds the available RF sequence.');
RF=sequence.RF(1:2:2*L); TR=TR_ms*ones(1,L);
zeta=mod(0:L-1,r)+1;
% Simulate a unit-density complex response for every pixel and time frame.
B=Blochmap_vectorized(RF,TR,truth(:,:,1),truth(:,:,2),0);
% Allocate full k-space storage; only measured rows receive signal and noise.
data=complex(zeros(n,m,L));
% Multiply by density before the unitary FFT. Simulate each time frame separately.
for t=1:L
    F=fft2(truth(:,:,3).*B(:,:,t))/sqrt(n*m);
    % Keep every r-th Fourier row. Independent real/imaginary Gaussian draws
    % give expected squared complex noise magnitude 2*noise_std^2 per sample.
    rows=zeta(t):r:n;
    data(rows,:,t)=F(rows,:)+noise_std*(randn(numel(rows),m)+1i*randn(numel(rows),m));
end
if use_discrepancy
    if isempty(noise_level)
        % sqrt(E||noise||^2) for independent real/imaginary N(0,sigma^2)
        % on measured entries only. This is an estimate, not a strict bound
        % on every Gaussian realization. Real data need a calibrated delta.
        noise_level=noise_std*sqrt(2*(n*m/r)*L);
    end
else
    noise_level=[]; %#ok<UNRCH> % Reached when the editable switch is false.
end

%% 3. Reference background initialization and projected LM reconstruction
reference_start=repmat(reshape(background_values,1,1,3),n,m,1);
% Fill tissue with uniform initial values while the background starts at
% reference T1/T2 and zero density. The mask imposes no constraint afterwards.
for j=1:3
    initial_map=reference_start(:,:,j); initial_map(mask)=initial_values(j);
    reference_start(:,:,j)=initial_map;
end
% Background pixels start at reference T1/T2 and zero density but can move.
options=struct('Scales',scales,'ReferenceTimes',reference_times,'H',1, ...
    'GaussNewton',gauss_newton,'MaxIterations',max_iterations, ...
    'Tolerance',stationarity_tolerance,'StepTolerance',step_tolerance, ...
    'BodyMask',mask,'LowerBounds',lower_bounds,'UpperBounds',upper_bounds, ...
    'InitialDamping',initial_damping,'DampingDecay',damping_decay, ...
    'ResidualDampingWeight',residual_damping_weight,'MinDamping',min_damping, ...
    'NoiseLevel',noise_level,'DiscrepancyFactor',discrepancy_factor, ...
    'LinearTolerance',linear_tolerance,'LinearMaxIterations',linear_max_iterations,'Verbose',true);
if ~isempty(parameters)
    % Copy all saved solver controls; retain the runtime maps/mask above.
    for field=fieldnames(parameters.lm_options).'
        name=field{1}; options.(name)=parameters.lm_options.(name);
    end
    options.H=parameters.h;
end
fprintf('Projected LM: %dx%d, L=%d, r=%d, physical L2 damping, data only, GN=%s\n', ...
    n,m,L,r,options.GaussNewton);
[reconstruction,info]=projected_lm_reconstruction(data,RF,TR,r,zeta,reference_start,options);
% The solver returns physical maps and histories that include the starting point.
% A budget, discrepancy or small-step stop is distinct from stationarity.
fprintf('Stopped by %s after %d updates / %d linear solves, %.2f seconds.\n', ...
    info.stop_reason,info.iterations,info.total_linear_solves,info.elapsed_seconds);
fprintf('Final outer projected gradient: %.6g; stationarity tolerance: %.6g.\n', ...
    info.projected_gradient(end),options.Tolerance);
if use_discrepancy
    fprintf('Measured residual: %.6g; tau*delta: %.6g; discrepancy satisfied: %d.\n', ...
        info.data_residual_norm(end),info.discrepancy_threshold,info.discrepancy_satisfied);
end
if strcmp(info.stop_reason,'linear_solve_failed')
    fprintf('Last linear flag=%d, relative residual=%.3g, lambda=%.3g.\n', ...
        info.last_linear_exitflag,info.last_linear_relative_residual,info.last_damping);
end

%% 4. Foreground errors and background deviations, using the physical maps
% Relative errors are dimensionless percentages; RMSE retains ms/density units.
% Background RMS here measures deviation from [reference_times,0].
all_results={reference_start,reconstruction}; labels={'Initial','Projected LM'};
relative_inside=zeros(2,3); rmse_inside=zeros(2,3); background_deviation_rms=zeros(2,3);
for k=1:2
    for j=1:3
        result_map=all_results{k}(:,:,j); target=truth(:,:,j);
        difference=result_map(mask)-target(mask);
        relative_inside(k,j)=100*norm(difference)/norm(target(mask));
        rmse_inside(k,j)=sqrt(mean(difference.^2));
        background_deviation_rms(k,j)=sqrt(mean((result_map(~mask)-background_values(j)).^2));
    end
end
fprintf('\nRelative error inside body (percent):\n');
disp(array2table(relative_inside,'VariableNames',{'T1','T2','Density'},'RowNames',labels));
fprintf('RMSE inside body (ms, ms, density units):\n');
disp(array2table(rmse_inside,'VariableNames',{'T1','T2','Density'},'RowNames',labels));
fprintf('RMS background deviations (T1-T1ref, T2-T2ref, density):\n');
disp(array2table(background_deviation_rms, ...
    'VariableNames',{'T1Deviation','T2Deviation','Density'},'RowNames',labels));

%% 5. Shared map/error reports, zoom insets, compact rows and profiles
method_names={'LM'};
export_reconstruction_report(truth,{reconstruction},reference_start,mask, ...
    method_names,design,output_folder);

%% 6. Convergence and solver diagnostics (independently callable plot routines)
fig=plot_lm_convergence(info,design);
export_reconstruction_pdf(fig,fullfile(output_folder,'convergence.pdf'),design);
fig=plot_lm_convergence(info,design,'diagnostics');
export_reconstruction_pdf(fig,fullfile(output_folder,'convergence.pdf'),design,true);
