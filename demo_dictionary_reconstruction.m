%% Standalone orthogonal dictionary-learning / denoising example
% Run demo_dictionary_reconstruction. This experiment uses noisy image patches (T1 by default);
% it does not require BLIP or the qMRI optimizer. Optional artificial maps are saved below.
% Workflow: corrupt a T1 image -> extract scaled patches -> learn D/C ->
% average reconstructed patches -> compare images and the stopping quantity.
% The noisy image is the training data; the clean image is used only for scoring.
% Reset the workspace/figures, then fix the random draw for reproducibility.
clearvars;
close all;
rng(1);
project_root = fileparts(mfilename('fullpath'));
% Build paths relative to this script so the helper folder is found reliably.
addpath(fullfile(project_root, 'helper_functions'));

%% Input: change the file and image variable here
% This denoising demo needs only one 2-D numeric image, not all three qMRI maps.
% For example, save('my_image.mat','T1'), then enter that file's path below.
input_file = fullfile(project_root, 'data', 'PWbrain256.mat');
% OPTIONAL artificial input: uncomment ALL THREE lines below to use it.
% Higher complexity adds more/smaller regions and finer texture. The MAT file is regenerated each run.
% [T1,T2,Density]=create_artificial_ground_truth(10,256); % Complexity, image size.
% input_file=fullfile(project_root,'data','artificial_maps.mat');
% save(input_file,'T1','T2','Density');
image_variable = 'T1';        % Name of the image array saved inside your MAT file.
% DOWNSAMPLING a 256-by-256 input: 1 -> 256-by-256, 2 -> 128-by-128, 4 -> 64-by-64.
pixel_stride = 1;             % Keep every pixel_stride-th row/column; no interpolation.
% Optional direct division is shown after loading the image below.
image_scale = 5000;           % Typical image amplitude, used only for numerical scaling.
noise_std = 150;              % Noise amplitude in the same units as your image.
% If you divide the image below, adjust image_scale and noise_std manually
% as appropriate. This denoising algorithm has no physical bounds.

%% 1. Learning settings and noisy image
patch_size = [8 8];           % All valid overlapping patches, without padding.
beta = 0.08;                 % L1 sparsity weight in the normalized patch units.
fixed_steps = 20;
tolerance = 0.1;             % Absolute combined change in D and C (paper's eta).
max_additional_steps = 1500;   % Cap for the tolerance-based continuation.

maps = load(input_file, image_variable);
truth = maps.(image_variable);
validateattributes(truth, {'numeric'}, {'real','finite','nonempty','2d'});
% Subsample by indexing, add real Gaussian image noise, then normalize
% amplitudes so beta and the D/C stopping tolerance act in dimensionless units.
truth = double(truth(1:pixel_stride:end, 1:pixel_stride:end));
% OPTIONAL: uncomment this line to divide the ground truth before adding noise.
% Adjust image_scale and noise_std above manually for your experiment.
% truth = truth / 10;
noisy = truth + noise_std*randn(size(truth));
normalized_image = noisy/image_scale;

%% Plot design: same defaults as demo_compare_all_reconstructions
% This demo has one image and no body mask; display and errors use all pixels.
display_lower=min(0,min(truth(:))); display_upper=max(truth(:));
if display_upper<=display_lower, display_upper=display_lower+1; end
design=reconstruction_plot_design(display_lower,display_upper,[1 3 4]);
design.MapNames={image_variable};
design.ShowFigures=true; design.ShowRowFigures=false;
design.ZeroBackground=true; design.RowShowTitles=false;
design.ProfileRow=[]; design.ZoomBox=[];
% MapLimits, ErrorLimits, MapColors, ErrorColors and sizes are editable in design.
output_folder=fullfile(project_root,'data','demo_reports',['dictionary_denoising_' datestr(now,'yyyymmdd_HHMMSS')]);

%% 2. Extract patches once; the dictionary function never sees ground truth
[R_mat, ~, overlap_count] = create_R_matrix(size(noisy), patch_size, 1);
% R_mat copies image pixels into all stride-one patches. X has one patch per
% column, so its size is prod(patch_size)-by-number_of_patches.
X = create_data_matrix_fast(normalized_image, patch_size, R_mat);
% Each pixel occurs in a different number of valid patches near the edges.
% These counts allow overlap averaging by elementwise division, with no
% dense diagonal matrix. Keep X and R_mat outside the inner learning loop.

%% 3. Fixed number of steps, followed by a warm-started accuracy-based solve
options = struct('StopRule','fixed', 'MaxIterations',fixed_steps, ...
    'LambdaD',1, 'LambdaC',1);
[D_fixed, C_fixed, fixed_info] = dictionary_reconstruction(X, beta, options);
% D_fixed is a square orthogonal patch basis; C_fixed holds sparse patch codes.

options.StopRule = 'tolerance';
% Continue from the computed D/C pair. Tolerance bounds the absolute joint
% D/C change per update; MaxIterations remains a cap if that test is not met.
options.MaxIterations = max_additional_steps;
options.Tolerance = tolerance;
options.D0 = D_fixed;
options.C0 = C_fixed;
[D, C, info] = dictionary_reconstruction(X, beta, options);

%% 4. Reconstruct images by averaging the learned patches
fixed_patches = D_fixed*C_fixed;
% R_mat' sums overlapping patch estimates back at their source pixels.
% Divide by coverage counts to average, reshape to an image, and undo scaling.
fixed_image = image_scale*reshape( ...
    (R_mat'*fixed_patches(:))./overlap_count, size(noisy));
learned_patches = D*C;
reconstruction = image_scale*reshape( ...
    (R_mat'*learned_patches(:))./overlap_count, size(noisy));

%% 5. Evaluate against ground truth, used only for these error measurements
images = {noisy, fixed_image, reconstruction};
relative_error_percent = zeros(3, 1);
rmse = zeros(3, 1);
for j = 1:3
    % Compare noisy, fixed-count and continued outputs with the same clean image.
    difference = images{j} - truth;
    relative_error_percent(j) = 100*norm(difference, 'fro')/norm(truth, 'fro');
    rmse(j) = sqrt(mean(difference(:).^2));
end
error_table = table(relative_error_percent, rmse, ...
    'VariableNames', {'RelativeError_percent','RMSE'}, ...
    'RowNames', {'Noisy','Fixed_steps','After_continuation'});
disp(error_table);
fprintf('Fixed solve: %d steps, %.3f s.\n', fixed_info.iterations, fixed_info.elapsed_seconds);
fprintf('Continuation: %d steps, %.3f s, stopped by %s.\n', ...
    info.iterations, info.elapsed_seconds, info.stop_reason);
fprintf('Final iterate change: %.4g; requested tolerance: %.4g.\n', ...
    info.step_norm(end), tolerance);
fprintf('Objective: %.6g -> %.6g. Nonzero coefficients: %.1f%%.\n', ...
    fixed_info.objective(1), info.objective(end), 100*info.nonzero_fraction(end));
% A smaller learning objective does not guarantee smaller ground-truth error.
% beta controls denoising versus detail loss. Ground truth never enters the
% learning updates or the stopping criterion.

%% 6. Shared image/error reports, zoom insets, compact rows and profiles
method_names={'Noisy','Fixed steps','After continuation'};
export_reconstruction_report(truth,images,noisy,true(size(truth)), ...
    method_names,design,output_folder);

%% 7. Dictionary stopping quantity
fig=reconstruction_report_figure('Dictionary convergence',design);
semilogy([fixed_info.step_norm; info.step_norm],'Color',design.MethodColors(3,:), ...
    'LineWidth',design.LineWidth);
% Join the histories; distinguish the continuation point from the tolerance.
hold on;
yline(tolerance,'--','Tolerance');
xline(fixed_info.iterations,':','Continuation');
xlabel('D/C update pairs'); ylabel('Combined iterate change'); grid on;
title('Paper stopping criterion');
export_reconstruction_pdf(fig,fullfile(output_folder,'convergence.pdf'),design);
