function parameters=reconstruction_parameters_128()
% RECONSTRUCTION_PARAMETERS_128 Saved qMRI comparison parameters.
% Usage: parameters=reconstruction_parameters_128();
% 128-by-128 baseline: settings from demo_compare_all_reconstructions.
% Edit this file to change the saved 128-by-128 settings.
% The two preset files are independent; neither loads the other.
% Sizes assume a 256-by-256 source. Input choice and plots stay in the demos.

% DOWNSAMPLING a 256-by-256 input: 1 -> 256-by-256, 2 -> 128-by-128, 4 -> 64-by-64.
pixel_stride=2;             % Keep every pixel_stride-th row/column; no interpolation.
% Optional map rescaling is selected in the demo's input-loading section.
% Bounds are physical [T1(ms), T2(ms), Density], chosen for your input maps.
% Keep time lower bounds positive and density lower bound zero.
lower_bounds=[100 10 0]; upper_bounds=[5500 600 110];
reference_times=[1000 100];  % Choose positive background times inside these bounds.
scales=[5000 500 100];       % Numerical scales in the SAME units; not intensity limits.
% Coordinate rule: x=(u-[reference_times,0])./scales. Changing scales changes
% the dictionary/spatial prior, so review Alpha/Lambda/Beta when changing units.
% If you rescale a map in the demo, adjust bounds, scales and reference_times here.
% Also review the BLIP grid spacing below for the new range.

%% Shared experiment and coarse-grid BLIP initialization
L=100; r=16; TR_ms=12;       % r must divide the reduced row count; r=1 works for any size.
noise_std=1.5; seed=1;       % Standard deviation of EACH real/imaginary component.
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

%% Algorithm settings (the demo adds shared coordinates/bounds at runtime)
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
% InitialDamping=r^2 follows Dong et al. (2019), sec. 4.4 (their s equals r).
% Decay and iteration budget are the settings checked for this implementation.
lm_options=struct('MaxIterations',15,'GaussNewton','approximate', ...
    'InitialDamping',r^2,'DampingDecay',.15,'MinDamping',1e-8, ...
    'ResidualDampingWeight',1e-15,'LinearTolerance',1e-8, ...
    'LinearMaxIterations',500,'Tolerance',1e-6,'StepTolerance',0, ...
    'NoiseLevel',[],'DiscrepancyFactor',1.1,'Verbose',true);
% NoiseLevel=[] disables discrepancy stopping. To use
% discrepancy stopping instead, set the noise norm after data generation.

% One-step DL: exactly one D/C update pair per map and outer iteration.
% Alpha weights spatial gradients; Lambda weights each map's patch fit;
% Lambda.*Beta weights coefficient sparsity. LambdaD/LambdaC control the
% proximal D/C updates. MapTolerance AND DictionaryTolerance must be met
% to stop for a small outer change. QP tolerances control each map subproblem.
one_options=struct( ...
    'MaxIterations',100,'Alpha',0.7, 'Lambda',[10.75 37.5 37.5]*1.5,'Beta',[.009 .009 .009]*1.5,...
    'PatchSize',[8 8],'StepNorm','H01','GaussNewton','approximate', ...
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

% Return plain MATLAB values; loading a preset does not simulate or save data.
parameters=struct();
parameters.pixel_stride=pixel_stride;
parameters.lower_bounds=lower_bounds;
parameters.upper_bounds=upper_bounds;
parameters.reference_times=reference_times;
parameters.scales=scales;
parameters.L=L;
parameters.r=r;
parameters.TR_ms=TR_ms;
parameters.noise_std=noise_std;
parameters.seed=seed;
parameters.h=h;
parameters.initial_use_body_mask=initial_use_body_mask;
parameters.coarse_spacing=coarse_spacing;
parameters.fine_spacing=fine_spacing;
parameters.coarse_T1=coarse_T1;
parameters.coarse_T2=coarse_T2;
parameters.initial_options=initial_options;
parameters.blip_options=blip_options;
parameters.lm_options=lm_options;
parameters.one_options=one_options;
parameters.nested_options=nested_options;
end
