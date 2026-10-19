%% Compare BLIP, projected LM, one-step DL and nested DL (MATLAB R2021b)
% Choose the input and saved parameter preset below, then Run. Requires Optimization Toolbox for DL.
% All four methods receive the SAME measurements and coarse-BLIP initial maps.
% T1/T2 use milliseconds; optional input rescaling changes the simulated phantom.
% Four reports, eight compact PDFs and results.mat are written to a new run folder in data.
% Workflow: define settings -> simulate one acquisition -> compute one common
% start -> run BLIP/LM/one-step DL/nested DL -> save metrics and PDF comparisons.
% truth is used for synthetic data and evaluation. The solvers receive data
% and initial maps; monitoring callbacks cannot change their update equations.
% Reset previous variables/figures before constructing this experiment.
clearvars;
close all;
root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'plot_routines'));
% Resolve the helper folders relative to this script, independent of the current folder.
addpath(fullfile(root,'helper_functions'),fullfile(root,'bloch_simulations'));

%% Input: choose reference maps here; physical ranges are in the parameter preset
% Supply a MAT file with equally sized 2-D arrays T1, T2 and Density.
% T1/T2 are in milliseconds; Density is nonnegative, with zero outside the body.
% To prepare your own file: save('my_maps.mat','T1','T2','Density').
% These maps are ground truth for SIMULATED measurements, not measured k-space.
%input_file=fullfile(root,'data','PWbrain256.mat'); % Or the full path to your MAT file.
% Artificial input (current default): comment ALL THREE lines to use PWbrain256 above.
% Higher complexity adds more/smaller regions and finer texture. The MAT file is regenerated each run.
[T1,T2,Density]=create_artificial_ground_truth(10,256); % Complexity, image size.
input_file=fullfile(root,'data','artificial_maps.mat');
save(input_file,'T1','T2','Density');

%% 1-2. Load saved experiment and algorithm parameters
% Choose ONE preset here. Edit its .m file to save changes for future runs.
% The 256 preset stores the standard full-resolution comparison settings.
% The two presets can be tuned independently.
parameters=reconstruction_parameters_128(); % Default: 128-by-128 from a 256 input.
% parameters=reconstruction_parameters_256(); % Uncomment to use full resolution.
% Optional one-run overrides go HERE, before the assignments below, e.g.:
% parameters.lm_options.MaxIterations=20;
% Bounds/scales/reference times and BLIP grids are configured in the preset.
% Changing map units requires matching changes to those settings and DL weights.
pixel_stride=parameters.pixel_stride;
lower_bounds=parameters.lower_bounds;
upper_bounds=parameters.upper_bounds;
reference_times=parameters.reference_times;
scales=parameters.scales;
L=parameters.L;
r=parameters.r;
TR_ms=parameters.TR_ms;
noise_std=parameters.noise_std;
seed=parameters.seed;
h=parameters.h;
initial_use_body_mask=parameters.initial_use_body_mask;
coarse_T1=parameters.coarse_T1;
coarse_T2=parameters.coarse_T2;
initial_options=parameters.initial_options;
blip_options=parameters.blip_options;
lm_options=parameters.lm_options;
one_options=parameters.one_options;
nested_options=parameters.nested_options;

%% 3. Plot design: shared color limits make the maps visually comparable
% These settings affect rendering, not reconstruction or numerical errors.
% Each row of MapLimits/ErrorLimits belongs to T1, T2 or density in that order.
design.MapColors={parula(256),parula(256),parula(256)}; % Replace by gray(256), etc.
design.MapLimits=[zeros(3,1) upper_bounds(:)];     % Display ranges follow physical bounds.
design.ErrorColors=hot(256);
design.ErrorLimits=[zeros(3,1) .2*(upper_bounds-lower_bounds).']; % Display first 20% of each range.
design.MethodColors=[.49 .18 .56; .47 .67 .19; .85 .33 .10; 0 .45 .74];
design.LineStyles={'-','--','-.','-'};             % BLIP, LM, one-step, nested.
design.FontName='Arial'; design.FontSize=10; design.LineWidth=1.4;
design.ReportFontName='Consolas'; design.ReportFontSize=10; % Monospaced summary table.
design.FigureSize=[1500 850]; design.Resolution=200;
design.ShowFigures=true; design.ZeroBackground=true; % Display copies only.
design.ProfileRow=[];                            % [] chooses the middle row.
design.ZoomBox=[]; % [row_first row_last col_first col_last] on reduced grid.
% [] uses the same automatic body region in maps.pdf and the eight compact PDFs.
% An explicit box also adds enlarged maps to profiles_and_summary.pdf.
design.ZoomInsetFraction=.43; % Maximum inset side relative to the full map.
design.ZoomColor=[1 .12 .28];  % ROI box, connector and inset border.
design.ComparisonFigureSize=[1100 1320]; % Portrait map/error comparison page.
design.RowPanelSize=300; design.RowPanelGap=14; % Compact PDF map size/gap in pixels.
design.RowShowTitles=false; % false removes headings for an image-only comparison row.
design.ShowRowFigures=false; % Export the eight compact PDFs without opening windows; true also requires ShowFigures=true.
output_folder=fullfile(root,'data',['comparison_all_' datestr(now,'yyyymmdd_HHMMSS')]); % PDF DESTINATION: a new time-stamped folder inside this project's data directory.
% The time-stamped folder holds all reports, compact rows and numerical results.

%% 4. Generate ONE data set in the current physical units
rng(seed);
% Reusing the seed and acquisition settings gives the same synthetic noise.
maps=load(input_file,'T1','T2','Density');
% Check the three arrays once and convert integer/single input to double arithmetic.
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
%warning: 'You are using a modified ground truth. Adjust bounds, scales and reference_times manually.'!!


outside=maps.Density==0;
% Positive background times keep the Bloch model valid. Their signal is still
% zero because density is zero there, so this does not introduce tissue signal.
maps.T1(outside)=reference_times(1); maps.T2(outside)=reference_times(2);
truth=cat(3,maps.T1(1:pixel_stride:end,1:pixel_stride:end), ...
    maps.T2(1:pixel_stride:end,1:pixel_stride:end), ...
    maps.Density(1:pixel_stride:end,1:pixel_stride:end));
mask=truth(:,:,3)>0; [n,m,~]=size(truth);
% Catch incompatible custom maps before expensive Bloch simulation/reconstruction.
assert(any(mask(:)),'The reduced input must contain positive Density values.');
for j=1:3
    values=truth(:,:,j);
    assert(all(values(mask)>=lower_bounds(j) & values(mask)<=upper_bounds(j)), ...
        'Input map %d is outside the physical bounds. Edit lower_bounds/upper_bounds.',j);
end
% assert stops with its message if a requested experiment is inconsistent.
% Require a valid sampling stride and compatible physical bounds/dictionary grids.
assert(r>=1 && r==fix(r) && mod(n,r)==0,'r must divide the image row count.');
assert(all(lower_bounds(1:2)>0) && all(upper_bounds>lower_bounds) && ...
    lower_bounds(3)==0,'Use positive time bounds and zero density lower bound.');
assert(all(reference_times>=lower_bounds(1:2) & reference_times<=upper_bounds(1:2)), ...
    'The reference times must lie inside the physical bounds.');
% Reject out-of-bounds grids instead of silently changing a requested experiment.
for choice={initial_options,blip_options}
    o=choice{1};
    assert(all(o.T1Grid>=lower_bounds(1) & o.T1Grid<=upper_bounds(1)) && ...
        all(o.T2Grid>=lower_bounds(2) & o.T2Grid<=upper_bounds(2)), ...
        'Both BLIP grids must lie inside the common physical bounds.');
end
sequence=load(fullfile(root,'data','sequence_parameter_middle.mat'),'RF');
% Select L pulses across the full stored sequence; TR is prescribed separately.
[RF,rf_indices]=select_rf_pulses(sequence.RF,L);
TR=TR_ms*ones(1,L);
zeta=mod(0:L-1,r)+1;
B=Blochmap_vectorized(RF,TR,truth(:,:,1),truth(:,:,2),0);
% B is n-by-m-by-L. Density multiplies the unit-density response before the FFT.
data=complex(zeros(n,m,L)); noise_squared=0;
for t=1:L
    % Simulate only acquired Fourier coefficients, with independent real and
    % imaginary Gaussian noise. Unmeasured array entries remain zero-filled.
    prediction=fft2(truth(:,:,3).*B(:,:,t))/sqrt(n*m);
    rows=zeta(t):r:n;
    noise=noise_std*(randn(numel(rows),m)+1i*randn(numel(rows),m));
    data(rows,:,t)=prediction(rows,:)+noise;
    noise_squared=noise_squared+sum(abs(noise(:)).^2);
    % This sum tracks the actual noise realization, not its expected norm.
end
noise_norm=sqrt(noise_squared); % Actual simulated noise; used only for diagnostic lines.
expected_noise_norm=noise_std*sqrt(2*n*m*L/r); % RMS norm inferred from noise_std.
% Optional early stopping examples (uncomment separately for either method):
% lm_options.NoiseLevel=expected_noise_norm;
% blip_options.NoiseLevel=expected_noise_norm;
clear B prediction noise maps;
% Release signal-simulation arrays; reconstruction still has data and truth for scoring.

%% 5. One coarse BLIP run, followed by the same initialization for everyone
initial_timer=tic;
u_zero=repmat(reshape([reference_times 0],1,1,3),n,m,1);
% Positive times and zero density give zero initial magnetization everywhere.
initial_options.DensityBounds=[lower_bounds(3) upper_bounds(3)];
fprintf('\nCommon coarse BLIP initialization: %dx%d, L=%d, r=%d\n',n,m,L,r);
[u0,initial_history]=blip_reconstruction(data,RF,TR,r,zeta,u_zero,initial_options);
% Save one common start; do not regenerate or refine it separately per method.
if initial_use_body_mask
    % Apply known support only to u0. Background pixels remain free in every solve.
    for j=1:3
        image=u0(:,:,j); zero_image=u_zero(:,:,j);
        image(~mask)=zero_image(~mask); u0(:,:,j)=image;
    end
end
initial_seconds=toc(initial_timer);
% This callback observes physical maps; truth never enters the update equations.
monitor=@(u) body_errors(u,truth,mask);
shared=struct('Scales',scales,'ReferenceTimes',reference_times,'H',h, ...
    'LowerBounds',lower_bounds,'UpperBounds',upper_bounds,'BodyMask',mask,'MonitorFcn',monitor);
blip_options.DensityBounds=[lower_bounds(3) upper_bounds(3)];
blip_options.MonitorFcn=monitor;
for field=fieldnames(shared).'
    % Dynamic field names copy the SAME scales, bounds and diagnostics into all options.
    name=field{1};
    lm_options.(name)=shared.(name);
    one_options.(name)=shared.(name); nested_options.(name)=shared.(name);
end
% Both dictionary solvers construct their PCA D0/C0 from exactly this u0.
% This deterministic setup is included in EACH solver's elapsed time.

%% 6. Run the four algorithms; keep every actual stopping reason
method_names={'BLIP','LM','One-step DL','Nested DL'};
method_options={blip_options,lm_options,one_options,nested_options};
reconstructions=cell(1,4); histories=cell(1,4); seconds=zeros(4,1);
% Cells hold map arrays and histories whose lengths may differ by stopping reason.
for method=1:4
    fprintf('\n--- %s ---\n',method_names{method}); method_timer=tic;
    if method==1
        [reconstructions{method},histories{method}]= ...
            blip_reconstruction(data,RF,TR,r,zeta,u0,method_options{method});
    elseif method==2
        [reconstructions{method},histories{method}]= ...
            projected_lm_reconstruction(data,RF,TR,r,zeta,u0,method_options{method});
    else
        % One-step and nested DL share the same outer solver; InnerStopRule differs.
        [reconstructions{method},~,~,histories{method}]= ...
            nested_dictionary_reconstruction(data,RF,TR,r,zeta,u0,method_options{method});
    end
    seconds(method)=toc(method_timer); state=histories{method};
    fprintf('Stopped by %s after %d steps, %.2f seconds.\n', ...
        state.stop_reason,state.iterations,seconds(method));
    if strcmp(state.stop_reason,'inner_tolerance_not_met')
        fprintf('Last inner changes [%g %g %g]; required eta=%g.\n', ...
            state.last_inner_step_norm,state.last_inner_tolerance);
    end
end
% Report separate solve times and totals including the shared coarse BLIP run.
% These are single-run timings, including setup/monitoring, excluding plotting.
relative_error=zeros(4,3); rmse=zeros(4,3); background_density=zeros(4,1);
% Final scalar metrics use the stored physical maps, before any display masking.
residual_norm=zeros(4,1); iterations=zeros(4,1); stop_reason=cell(4,1);
for method=1:4
    u=reconstructions{method}; state=histories{method};
    relative_error(method,:)=state.monitor(end,:);
    for j=1:3
        difference=u(:,:,j)-truth(:,:,j); rmse(method,j)=sqrt(mean(difference(mask).^2));
    end
    density=u(:,:,3); background_density(method)=sqrt(mean(density(~mask).^2));
    % DL stores h^2*||residual||^2/2, whereas BLIP/LM store the residual norm.
    % Convert them to the same measured-data quantity for fair comparison.
    if method<=2, residual=state.data_residual_norm; else, residual=sqrt(2*state.data)/h; end
    histories{method}.comparison_residual=residual;
    residual_norm(method)=residual(end); iterations(method)=state.iterations;
    stop_reason{method}=state.stop_reason;
end
total_seconds=seconds+initial_seconds;
% Each total includes the shared coarse initialization once plus that method's solve.
comparison=table(iterations,seconds,total_seconds,relative_error(:,1),relative_error(:,2), ...
    relative_error(:,3),residual_norm,background_density,stop_reason, ...
    'VariableNames',{'Steps','SolveSeconds','TotalSeconds','T1Error_percent','T2Error_percent', ...
    'DensityError_percent','DataResidual','BackgroundDensityRMS','StopReason'},'RowNames',method_names);
disp(comparison);
mkdir(output_folder);
% Only now create the output directory, after the numerical solves have completed.
% Save enough to inspect the settings/results later without rerunning the solvers.
experiment=struct('InputFile',input_file,'PixelStride',pixel_stride,'Frames',L,'Undersampling',r,'TR_ms',TR_ms, ...
    'RFSelection','uniform_full_sequence','RFIndices',rf_indices,'RFSourceLength',numel(sequence.RF), ...
    'NoiseStd',noise_std,'Seed',seed,'Scales',scales,'ReferenceTimes',reference_times, ...
    'LowerBounds',lower_bounds,'UpperBounds',upper_bounds,'H',h, ...
    'InitialUseBodyMask',initial_use_body_mask,'InitialSeconds',initial_seconds, ...
    'NoiseNorm',noise_norm,'ExpectedNoiseNorm',expected_noise_norm);
% Do not serialize local function handles; diagnostics themselves are already saved.
saved_options=method_options; saved_histories=histories;
% Work on copies so removing function handles for serialization leaves live histories intact.
for method=1:4
    saved_options{method}.MonitorFcn=[]; saved_histories{method}.options.MonitorFcn=[];
end
save(fullfile(output_folder,'results.mat'),'truth','mask','u0','reconstructions', ...
    'saved_histories','saved_options','initial_options','initial_history','comparison', ...
    'rmse','experiment','design','method_names','RF','TR','zeta');

%% 7. Plot/export the saved results; each figure is also callable separately
% Loading the saved structure exercises the same interface used for later replotting.
% Examples: plot_comparison(results), plot_convergence(results,'map-change').
export_comparison_reports(fullfile(output_folder,'results.mat'),output_folder);
fprintf('\nPDFs and results saved in:\n%s\n',output_folder);

%% Evaluation helper used by solver monitoring

function errors=body_errors(u,truth,mask)
    % Relative L2 errors on foreground pixels only; eps guards a zero reference norm.
    errors=zeros(1,3);
    for j=1:3
        target=truth(:,:,j); difference=u(:,:,j)-target;
        errors(j)=100*norm(difference(mask))/max(norm(target(mask)),eps);
    end
end
