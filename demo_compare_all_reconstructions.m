%% Compare BLIP, projected LM, one-step DL and nested DL (MATLAB R2021b)
% Edit the input section and sections 1-3, then Run. Requires Optimization Toolbox for DL.
% All four methods receive the SAME measurements and coarse-BLIP initial maps.
% T1/T2 use milliseconds; optional input rescaling changes the simulated phantom.
% Four PDFs and a small results MAT file are written to a new run folder in data.
% Workflow: define settings -> simulate one acquisition -> compute one common
% start -> run BLIP/LM/one-step DL/nested DL -> save metrics and PDF comparisons.
% truth is used for synthetic data and evaluation. The solvers receive data
% and initial maps; monitoring callbacks cannot change their update equations.
% Reset previous variables/figures before constructing this experiment.
clearvars;
close all;
root=fileparts(mfilename('fullpath'));
% Resolve the helper folders relative to this script, independent of the current folder.
addpath(fullfile(root,'helper_functions'),fullfile(root,'bloch_simulations'));

%% Input: change the file and physical ranges here
% Supply a MAT file with equally sized 2-D arrays T1, T2 and Density.
% T1/T2 are in milliseconds; Density is nonnegative, with zero outside the body.
% To prepare your own file: save('my_maps.mat','T1','T2','Density').
% These maps are ground truth for SIMULATED measurements, not measured k-space.
input_file=fullfile(root,'data','PWbrain256.mat'); % Or the full path to your MAT file.
% OPTIONAL artificial input: uncomment ALL THREE lines below to use it.
% Higher complexity adds more/smaller regions and finer texture. The MAT file is regenerated each run.
% [T1,T2,Density]=create_artificial_ground_truth(10,256); % Complexity, image size.
% input_file=fullfile(root,'data','artificial_maps.mat');
% save(input_file,'T1','T2','Density');
% DOWNSAMPLING a 256-by-256 input: 1 -> 256-by-256, 2 -> 128-by-128, 4 -> 64-by-64.
pixel_stride=2;             % Keep every pixel_stride-th row/column; no interpolation.
% Optional direct division is shown after loading the maps below.
% Bounds are physical [T1(ms), T2(ms), Density], chosen for your input maps.
% Keep time lower bounds positive and density lower bound zero.
lower_bounds=[100 10 0]; upper_bounds=[5500 600 110];
reference_times=[1000 100];  % Choose positive background times inside these bounds.
scales=[5000 500 100];       % Numerical scales in the SAME units; not intensity limits.
% Coordinate rule: x=(u-[reference_times,0])./scales. Changing scales changes
% the dictionary/spatial prior, so review Alpha/Lambda/Beta when changing units.
% If you divide a map below, adjust its bounds, scales and reference_times manually.
% Also review the BLIP grid spacing below for the new range.

%% 1. Shared experiment and coarse-grid BLIP initialization
L=100; r=16; TR_ms=12;       % r must divide the reduced row count; r=1 works for any size.
noise_std=1.0; seed=1;       % Standard deviation of EACH real/imaginary component.
h=1;
initial_use_body_mask=true; % Legacy initialization: zero density outside body.
                           % The mask is NOT a constraint in any subsequent solve.
coarse_spacing=[200 50];    % [T1 T2] grid spacing in ms; smaller costs more memory/time.
fine_spacing=[100 10];      % The reconstruction dictionary uses this finer spacing.
% Both grids follow the bounds above, including the upper endpoints and reference times.
coarse_T1=unique([lower_bounds(1):coarse_spacing(1):upper_bounds(1) upper_bounds(1) reference_times(1)]);
coarse_T2=unique([lower_bounds(2):coarse_spacing(2):upper_bounds(2) upper_bounds(2) reference_times(2)]);
initial_options=struct('T1Grid',coarse_T1,'T2Grid',coarse_T2, ...
    'MaxIterations',20,'MaxTrials',20,'InitialStep',r,'Kappa',.99, ...
    'MaxBacktracking',30,'BacktrackingFactor',2,'BatchSize',512, ...
    'StepTolerance',0,'Verbose',false);
% The legacy initializer spends 20 TRIALS, including rejected trials.
% Its grid is now in physical milliseconds and is freely adjustable above.

%% 2. Separate algorithm settings (shared coordinates/bounds are added below)
% BLIP: a finer Bloch dictionary. Including the coarse grid makes the common
% nonzero initial signals members of this dictionary as well.
blip_options=struct('T1Grid',unique([lower_bounds(1):fine_spacing(1):upper_bounds(1) coarse_T1]), ...
    'T2Grid',unique([lower_bounds(2):fine_spacing(2):upper_bounds(2) coarse_T2]), ...
    'MaxIterations',30,'MaxTrials',Inf,'InitialStep',r,'Kappa',.99, ...
    'MaxBacktracking',30,'BacktrackingFactor',2,'BatchSize',512, ...
    'StepTolerance',0,'NoiseLevel',[],'DiscrepancyFactor',1.1,'Verbose',true);

% LM: physical L2 damping, no spatial or dictionary regularization.
% GaussNewton chooses local Gram curvature or exact sampled FFT products.
% Tolerance tests map stationarity; LinearTolerance tests the normal equations.
% InitialDamping/DampingDecay/MinDamping define the geometric damping schedule.
lm_options=struct('MaxIterations',15,'GaussNewton','approximate', ...
    'InitialDamping',r^2,'DampingDecay',.15,'MinDamping',1e-8, ...
    'ResidualDampingWeight',1e-15,'LinearTolerance',1e-8, ...
    'LinearMaxIterations',500,'Tolerance',1e-6,'StepTolerance',0, ...
    'NoiseLevel',[],'DiscrepancyFactor',1.1,'Verbose',true);
% NoiseLevel=[] gives the paper's manually chosen iteration budget. To use
% discrepancy stopping instead, set the noise norm after data generation.

% One-step DL: exactly one D/C update pair per map and outer iteration.
% Alpha weights spatial gradients; Lambda weights each map's patch fit;
% Lambda.*Beta weights coefficient sparsity. LambdaD/LambdaC control the
% proximal D/C updates. MapTolerance AND DictionaryTolerance must be met
% to stop for a small outer change. QP tolerances control each map subproblem.
one_options=struct( ...
    'MaxIterations',100,'Alpha',0.5, 'Lambda',[6.75 37.5 37.5]*2.0,'Beta',[.009 .009 .009]*1.5,...
    'PatchSize',[6 6],'StepNorm','H01','GaussNewton','approximate', ...
    'InnerStopRule','one-step','LambdaD',1,'LambdaC',1, ...
    'MapTolerance',1e-6,'DictionaryTolerance',1e-4, ...
    'InitialDamping',1,'MinDamping',1e-9,'MaxBacktracking',20, ...
    'QPOptimalityTolerance',1e-7,'QPFunctionTolerance',0, ...
    'QPPCGTolerance',1e-6,'QPRequireOptimality',true,'QPMaxIterations',600, ...
    'Verbose', true);

% Nested DL: start with the same settings; override any field independently.
% Same PatchSize/Beta gives the same deterministic D0/C0 from the common u0.
nested_options=one_options;
% Copying the structure preserves all shared choices; following assignments
% override only the amount and accuracy of dictionary learning per outer step.
nested_options.InnerStopRule='tolerance';
nested_options.InnerMaxIterations=1000;
nested_options.InnerToleranceScale=.005;
nested_options.InnerDecay=.75; % eta_k=scale*||(D0,C0)||/(k+1)^gamma; gamma>1/2.
% nested_options.Lambda=[4.5 45 45]*3.5; % Example of a separate method setting.

%% 3. Plot design: shared color limits make the maps visually comparable
% These settings affect rendering, not reconstruction or numerical errors.
% Each row of MapLimits/ErrorLimits belongs to T1, T2 or density in that order.
design.MapColors={turbo(256),turbo(256),turbo(256)}; % Replace by gray(256), etc.
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
% [] uses an automatic body region for the insets on page 2 of maps.pdf.
% An explicit box also adds enlarged maps to profiles_and_summary.pdf.
design.ZoomInsetFraction=.43; % Maximum inset side relative to the full map.
design.ZoomColor=[1 .12 .28];  % ROI box, connector and inset border.
design.ComparisonFigureSize=[1100 1320]; % Portrait map/error comparison page.
output_folder=fullfile(root,'data',['comparison_all_' datestr(now,'yyyymmdd_HHMMSS')]);
% The time-stamped folder holds this run's four reports and numerical results.

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
% Select L stored RF pulses, prescribe TR, and cycle the sampling offset 1:r.
RF=sequence.RF(1:2:2*L); TR=TR_ms*ones(1,L); zeta=mod(0:L-1,r)+1;
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
    'rmse','experiment','design','method_names');

%% 7. Map and absolute-error PDFs (equal color limits within each row)
% First page: truth, common start and four results. Second page: linked zooms.
% All displays operate on copies, and color saturation does not clip stored values.
map_names={'T1 (ms)','T2 (ms)','Density'};
fig=new_figure('Parameter maps',design); tiledlayout(3,6,'TileSpacing','compact','Padding','compact');
shown=[{truth,u0} reconstructions]; headings=[{'Truth','Coarse BLIP start'} method_names];
for j=1:3
    for column=1:6
        ax=nexttile;
        show_map(ax,shown{column}(:,:,j),mask,design.MapLimits(j,:),design.MapColors{j},design);
        title(sprintf('%s\n%s',headings{column},map_names{j}));
        if column==6, colorbar; end
    end
end
sgtitle(sprintf('Shared start | %dx%d, L=%d, r=%d, noise SD=%.2g per component',n,m,L,r,noise_std));
export_pdf(fig,fullfile(output_folder,'maps.pdf'),design,false);
% Publication-style grid: truth and all methods, alternating map/error rows.
export_reconstruction_comparison(truth,reconstructions,mask,method_names, ...
    design,fullfile(output_folder,'maps.pdf'),true);

fig=new_figure('Absolute reconstruction errors',design);
% Absolute errors are in ms for T1/T2 and the density map's original units.
tiledlayout(3,4,'TileSpacing','compact','Padding','compact');
for j=1:3
    for method=1:4
        ax=nexttile;
        show_map(ax,abs(reconstructions{method}(:,:,j)-truth(:,:,j)),mask, ...
            design.ErrorLimits(j,:),design.ErrorColors,design);
        title(sprintf('%s | %s\nBody error %.2f%%',method_names{method},map_names{j},relative_error(method,j)));
        if method==4, colorbar; end
    end
end
sgtitle('Absolute errors in physical units; common limits per row (larger errors saturate)');
export_pdf(fig,fullfile(output_folder,'errors.pdf'),design,false);

%% 8. Convergence PDF: residuals, reconstruction accuracy, dictionary work
% Append three pages to one PDF: fit/runtime, body errors, and DL diagnostics.
% Curves of states start at zero; curves of accepted-step work start at one.
convergence_pdf=fullfile(output_folder,'convergence.pdf');
fig=new_figure('Data fit and runtime',design); tiledlayout(2,2,'TileSpacing','compact');
for panel=1:2
    % Use iteration count in panel 1 and elapsed solver time in panel 2.
    nexttile; hold on;
    for method=1:4
        state=histories{method}; horizontal=0:state.iterations;
        if panel==2, horizontal=state.time; end
        draw_curve(horizontal,state.comparison_residual,method,design);
    end
    if noise_norm>0, yline(noise_norm,':k','Actual noise norm','HandleVisibility','off'); end
    set(gca,'YScale','log'); grid on; ylabel('Measured residual norm');
    if panel==1, xlabel('Completed step'); else, xlabel('Solver seconds (including setup)'); end
    title('Same data-fit quantity for all methods'); legend(method_names,'Location','best');
end
nexttile; bar([seconds initial_seconds*ones(4,1)],'stacked');
% The stacked bars distinguish work unique to a method from the common start.
set(gca,'XTick',1:4,'XTickLabel',method_names); ylabel('Seconds'); grid on;
legend({'Solver + its setup','Shared coarse initialization'},'Location','best'); title('Single-run time costs');
nexttile; bar(background_density); set(gca,'XTick',1:4,'XTickLabel',method_names);
ylabel('RMS density outside body'); grid on; title('Background leakage (unmasked reconstruction)');
sgtitle('Data fit and computational cost; iteration costs differ across methods');
export_pdf(fig,convergence_pdf,design,false);

fig=new_figure('Errors over iterations and time',design); tiledlayout(2,3,'TileSpacing','compact');
% MonitorFcn stored three body errors at every state; plot them without rerunning solvers.
for row=1:2
    for j=1:3
        nexttile; hold on;
        for method=1:4
            state=histories{method}; horizontal=0:state.iterations;
            if row==2, horizontal=state.time; end
            draw_curve(horizontal,state.monitor(:,j),method,design);
        end
        grid on; title(map_names{j}); ylabel('Relative body error (%)');
        if row==1, xlabel('Completed step'); else, xlabel('Solver seconds (including setup)'); end
        if j==1, legend(method_names,'Location','best'); end
    end
end
sgtitle('Reconstruction quality, including the common initial maps at step zero');
export_pdf(fig,convergence_pdf,design,true);

fig=new_figure('Dictionary convergence',design); tiledlayout(2,2,'TileSpacing','compact');
% Full objectives are comparable between the DL methods only for equal prior weights.
% The map projected gradient treats D/C as fixed and does not certify joint optimality.
for panel=1:4
    nexttile; hold on;
    for method=3:4
        state=histories{method}; horizontal=0:state.iterations;
        switch panel
            case 1, vertical=state.objective; label='Full regularized objective';
            case 2, vertical=state.projected_gradient; label='Map projected gradient (D,C fixed)';
            case 3
                horizontal=1:state.iterations; vertical=sum(state.inner_iterations,2);
                label='D/C pairs per outer step (sum over maps)';
            case 4
                horizontal=1:state.iterations; vertical=state.qp_firstorderopt;
                label='QP first-order optimality';
        end
        draw_curve(horizontal,max(vertical,eps),method,design);
    end
    set(gca,'YScale','log'); grid on; xlabel('Completed step'); title(label);
    legend(method_names(3:4),'Location','best');
end
sgtitle('DL diagnostics: objective values are comparable only for equal regularization settings');
export_pdf(fig,convergence_pdf,design,true);

%% 9. Profiles and a compact numerical/settings summary
% Profiles compare one image row. An explicit ZoomBox also adds a crop page;
% the final page records settings and actual stopping reasons alongside errors.
profiles_pdf=fullfile(output_folder,'profiles_and_summary.pdf');
if isempty(design.ProfileRow), design.ProfileRow=round(n/2); end
row=design.ProfileRow;
fig=new_figure('Parameter profiles',design); tiledlayout(1,3,'TileSpacing','compact');
for j=1:3
    nexttile; hold on; target=truth(:,:,j);
    if design.ZeroBackground, target(~mask)=0; end
    plot(1:m,target(row,:),'k:','LineWidth',2);
    for method=1:4
        values=reconstructions{method}(:,:,j);
        if design.ZeroBackground, values(~mask)=0; end
        draw_curve(1:m,values(row,:),method,design);
    end
    grid on; xlabel('Column'); ylabel(map_names{j}); title(sprintf('Row %d',row));
    if j==1, legend([{'Truth'} method_names],'Location','best'); end
end
sgtitle('Physical parameter profiles (display background follows design.ZeroBackground)');
export_pdf(fig,profiles_pdf,design,false);
if ~isempty(design.ZoomBox)
    % Rows/columns refer to the reduced reconstruction grid and include both endpoints.
    box=design.ZoomBox; rows=box(1):box(2); columns=box(3):box(4);
    fig=new_figure('Enlarged parameter maps',design); tiledlayout(3,5,'TileSpacing','compact');
    shown=[{truth} reconstructions]; headings=[{'Truth'} method_names];
    for j=1:3
        for column=1:5
            ax=nexttile;
            show_map(ax,shown{column}(rows,columns,j),mask(rows,columns), ...
                design.MapLimits(j,:),design.MapColors{j},design);
            title([headings{column} ' | ' map_names{j}]); if column==5, colorbar; end
        end
    end
    sgtitle(sprintf('Zoom: rows %d-%d, columns %d-%d',box));
    export_pdf(fig,profiles_pdf,design,true);
end
fig=new_figure('Experiment summary',design); ax=axes(fig,'Position',[.04 .04 .92 .9]); axis(ax,'off');
lines_out={sprintf('FOUR-METHOD COMPARISON | %dx%d | L=%d | r=%d | noise SD=%.3g | seed=%d',n,m,L,r,noise_std,seed), ...
    sprintf('TR=%.3g ms; scales=[%g %g %g]; references=[%g %g]; h=%g',TR_ms,scales,reference_times,h), ...
    sprintf('Common coarse BLIP: %d x %d candidates; %d accepted / %d trials; %.2f s; mask applied=%d', ...
    numel(coarse_T1),numel(coarse_T2),initial_history.iterations,initial_history.total_trials,initial_seconds,initial_use_body_mask), ...
    ' ', 'Method       Steps  Solve(s)  Total(s)    T1(%)   T2(%)  Density(%)   Residual   Background rho RMS'};
for method=1:4
    lines_out{end+1}=sprintf('%-12s %4d %9.2f %9.2f %8.2f %7.2f %11.2f %10.3g %12.3g', ...
        method_names{method},iterations(method),seconds(method),total_seconds(method), ...
        relative_error(method,:),residual_norm(method),background_density(method));
end
lines_out=[lines_out {' ','Stopping reasons (a budget or failed inner solve does not certify convergence):'}];
for method=1:4
    lines_out{end+1}=sprintf('%-12s %s',method_names{method},stop_reason{method});
end
lines_out=[lines_out {' ',sprintf('BLIP: %d x %d candidates, max steps=%d; LM: %s GN, max steps=%d, damping=%g * %g^k', ...
    numel(blip_options.T1Grid),numel(blip_options.T2Grid),blip_options.MaxIterations, ...
    lm_options.GaussNewton,lm_options.MaxIterations,lm_options.InitialDamping,lm_options.DampingDecay)}];
for method=3:4
    o=method_options{method};
    lines_out{end+1}=sprintf('%s: alpha=%g; lambda=[%g %g %g]; beta=[%g %g %g]; patch=%dx%d; max steps=%d', ...
        method_names{method},o.Alpha,o.Lambda,o.Beta,o.PatchSize,o.MaxIterations);
    lines_out{end+1}=sprintf('  %s GN; %s step norm; QP tolerances (opt/function/PCG)=[%g %g %g]', ...
        o.GaussNewton,o.StepNorm,o.QPOptimalityTolerance,o.QPFunctionTolerance,o.QPPCGTolerance);
end
lines_out=[lines_out {sprintf('Nested eta: %.4g * ||(D0,C0)|| / (k+1)^%.3g; inner cap=%d', ...
    nested_options.InnerToleranceScale,nested_options.InnerDecay,nested_options.InnerMaxIterations), ...
    ' ', 'Errors use body pixels; background density RMS uses the unmasked output.', ...
    'Map/profile background masking is for display only. Full settings and RMSE are in results.mat.', ...
    'Full DL objectives include regularizers; BLIP/LM use data fit only. Compare their residuals above.', ...
    'Times are one run each and include monitoring/setup. Common initialization is added in Total(s).'}];
% Interpreter none preserves literal underscores and option names in the report.
text(ax,0,1,lines_out,'VerticalAlignment','top','Interpreter','none', ...
    'FontName',design.ReportFontName,'FontSize',design.ReportFontSize);
export_pdf(fig,profiles_pdf,design,true);
fprintf('\nPDFs and results saved in:\n%s\n',output_folder);

%% Small plotting/evaluation helpers; none changes the reconstructed maps
function errors=body_errors(u,truth,mask)
    % Relative L2 errors on foreground pixels only; eps guards a zero reference norm.
    errors=zeros(1,3);
    for j=1:3
        target=truth(:,:,j); difference=u(:,:,j)-target;
        errors(j)=100*norm(difference(mask))/max(norm(target(mask)),eps);
    end
end

function fig=new_figure(name,design)
    % Apply the same geometry/fonts to each report page; off still allows PDF export.
    visible='off'; if design.ShowFigures, visible='on'; end
    fig=figure('Name',name,'Color','w','Visible',visible, ...
        'Position',[40 40 design.FigureSize],'DefaultAxesFontName',design.FontName, ...
        'DefaultAxesFontSize',design.FontSize,'DefaultTextFontName',design.FontName, ...
        'DefaultTextFontSize',design.FontSize);
end

function show_map(ax,values,mask,limits,colors,design)
    % values is a plotting copy; masking and intensity limits leave solver maps untouched.
    if design.ZeroBackground, values(~mask)=0; end
    imagesc(ax,values,limits); axis(ax,'image'); axis(ax,'off'); colormap(ax,colors);
end

function draw_curve(x,y,method,design)
    % Reuse each method's color/style across every convergence and profile panel.
    plot(x,y,'Color',design.MethodColors(method,:),'LineStyle',design.LineStyles{method}, ...
        'LineWidth',design.LineWidth);
end

function export_pdf(fig,path,design,append)
    % Hide interactive controls, finish pending graphics updates, then write a page.
    axes_list=findall(fig,'Type','axes');
    for ax=axes_list.'
        if ~isempty(ax.Toolbar), ax.Toolbar.Visible='off'; end % Exclude hover controls.
    end
    drawnow;
    % append=false starts the PDF; append=true adds the next page to it.
    exportgraphics(fig,path,'ContentType','vector','BackgroundColor','white', ...
        'Resolution',design.Resolution,'Append',append);
    if ~design.ShowFigures, close(fig); end
end
