%% Projected Levenberg-Marquardt qMRI reconstruction (MATLAB R2021b)
% Uses the Bloch signal model followed by sampled unitary Fourier transforms.
% DATA-ONLY objective. L2 damping penalizes PHYSICAL increments, as in the
% pasted algorithm. Each update is an unconstrained solve followed by clipping.
% Damping follows a schedule; objective increases are retained.
% Known body support is used for initialization, diagnostics and display.
% Workflow: generate one noisy data set -> initialize physical maps -> solve
% damped normal equations and clip to bounds -> inspect fit, errors and steps.
% Clear previous variables/figures, then seed random noise for reproducible runs.
clearvars;
close all;
rng(1);
% Locate data and helpers relative to this script's own folder.
root=fileparts(mfilename('fullpath'));
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
L=80;
r=8;                      % Must divide the reduced row count; r=1 works for any size.
noise_std=2.2;             % Standard deviation of EACH real/imaginary part.
TR_ms=12;
max_iterations=15;
gauss_newton='approximate'; % Sparse direct solve; 'exact' uses sampled FFT + PCG.
initial_damping=r^2;      % Damping before the first update (k=0).
damping_decay=.15;        % lambda_k=max(floor,lambda_0*decay^k,weight*residual).
residual_damping_weight=1e-15;
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
fprintf('Projected LM: %dx%d, L=%d, r=%d, physical L2 damping, data only, GN=%s\n', ...
    n,m,L,r,gauss_newton);
[reconstruction,info]=projected_lm_reconstruction(data,RF,TR,r,zeta,reference_start,options);
% The solver returns physical maps and histories that include the starting point.
% A budget, discrepancy or small-step stop is distinct from stationarity.
fprintf('Stopped by %s after %d updates / %d linear solves, %.2f seconds.\n', ...
    info.stop_reason,info.iterations,info.total_linear_solves,info.elapsed_seconds);
fprintf('Final outer projected gradient: %.6g; stationarity tolerance: %.6g.\n', ...
    info.projected_gradient(end),stationarity_tolerance);
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

%% 5. Physical maps and a separate zero-background display
% Make both versions from copies so hiding the background cannot change metrics.
plot_names={'T1 (ms)','T2 (ms)','Density'};
for masked_display=[false true]
    figure_name='Projected LM: physical maps'; caption='Physical maps';
    if masked_display
        figure_name='Projected LM: zero-background display';
        caption='Background set to zero for display';
    end
    figure('Name',figure_name,'Color','w','Position',[50 50 1000 850]);
    tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
    for j=1:3
        target=truth(:,:,j); display_map=reconstruction(:,:,j);
        limits=[0 max(target(:))];
        if masked_display, target(~mask)=0; display_map(~mask)=0; end
        ax=nexttile; imagesc(target,limits); axis image off; colormap(ax,gray); colorbar;
        title([plot_names{j},': truth']);
        ax=nexttile; imagesc(display_map,limits); axis image off; colormap(ax,gray); colorbar;
        title(sprintf('%s: projected LM\nError in body: %.1f%%',plot_names{j},relative_inside(2,j)));
    end
    sgtitle(sprintf('%s; physical L2 damping; GN: %s',caption,gauss_newton));
end

%% 6. Data residual, stationarity, background density and unmasked profiles
% Include iteration zero. realmin avoids log(0) in plots without changing histories.
figure('Name','Projected LM: convergence and profiles','Color','w','Position',[80 80 1200 750]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
nexttile; semilogy(0:info.iterations,max(info.data_residual_norm,realmin),'LineWidth',1.2);
if use_discrepancy && info.discrepancy_threshold>0
    hold on; yline(info.discrepancy_threshold,'k--','tau*delta');
end
grid on; xlabel('Update'); ylabel('Residual norm'); title('Measured data and noise threshold');
nexttile; semilogy(0:info.iterations,max(info.projected_gradient,realmin),'LineWidth',1.2);
grid on; xlabel('Update'); ylabel('Normalized coordinates'); title('Projected gradient norm');
nexttile; plot(0:info.iterations,info.background_rms(:,3),'LineWidth',1.2);
grid on; xlabel('Update'); ylabel('RMS density'); title('Density leaking into background');
profile_row=round(n/2);
for j=1:3
    nexttile; plot(truth(profile_row,:,j),'k--','LineWidth',1.5); hold on;
    plot(reconstruction(profile_row,:,j),'LineWidth',1.1);
    grid on; xlabel('Column'); ylabel(plot_names{j}); title(sprintf('Row %d profile',profile_row));
end
legend({'Truth','Projected LM'},'Location','best');
sgtitle(sprintf('Projected LM; physical L2 damping; GN: %s',gauss_newton));

%% 7. Scheduled damping, linear accuracy and the effect of projection
% These histories have one value per completed update. Linear residuals are
% measured before clipping; step norms use the actual increment after clipping.
figure('Name','Projected LM: solver diagnostics','Color','w','Position',[100 100 1100 650]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile; semilogy(1:info.iterations,info.damping,'LineWidth',1.2);
grid on; xlabel('Update'); ylabel('Lambda'); title('Scheduled damping');
nexttile; semilogy(1:info.iterations,max(info.linear_relative_residual,realmin),'LineWidth',1.2);
hold on; yline(linear_tolerance,'k--'); grid on;
xlabel('Update'); ylabel('Relative residual'); title('Linear system before projection');
nexttile; plot(1:info.iterations,100*info.projection_fraction,'LineWidth',1.2);
grid on; xlabel('Update'); ylabel('Percent of entries'); title('Entries clipped to bounds');
nexttile; semilogy(1:info.iterations,max(info.step_norm,realmin),'LineWidth',1.2);
grid on; xlabel('Update'); ylabel('Relative step'); title('Actual step after projection');
