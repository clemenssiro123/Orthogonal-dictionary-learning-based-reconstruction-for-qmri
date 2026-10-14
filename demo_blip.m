%% Simple BLIP reconstruction experiment (MATLAB R2021b)
% Run this script from the Editor or type demo_blip in the Command Window.
% It creates noisy, undersampled Fourier data, calls the existing BLIP
% routine, and compares its output with the known parameter maps.
% Default inputs are read only; the optional artificial input is saved explicitly below.
% Workflow: choose acquisition settings -> load truth -> simulate sampled
% measurements -> build a Bloch dictionary -> reconstruct -> evaluate/plot.
% This demo calls the legacy helper BLIP.m; the four-method demo uses the
% newer blip_reconstruction interface with explicit histories and bounds.
% clearvars removes results of earlier runs; close all removes their figures.
clearvars;
close all;

project_root = fileparts(mfilename('fullpath'));
% Resolve helper/data folders relative to this script, not the current folder.
addpath(fullfile(project_root, 'bloch_simulations'));
addpath(fullfile(project_root, 'helper_functions'));

%% Input: change the file and relaxation-time ranges here
% MAT file with equally sized 2-D arrays T1, T2 (ms) and nonnegative Density.
% Density=0 marks background. Save your maps with save('my_maps.mat','T1','T2','Density').
% These are reference maps for simulated measurements; optional rescaling is selected below.
input_file = fullfile(project_root, 'data', 'PWbrain256.mat'); % Or your full MAT-file path.
% OPTIONAL artificial input: uncomment ALL THREE lines below to use it.
% Higher complexity adds more/smaller regions and finer texture. The MAT file is regenerated each run.
% [T1,T2,Density]=create_artificial_ground_truth(10,256); % Complexity, image size.
% input_file=fullfile(project_root,'data','artificial_maps.mat');
% save(input_file,'T1','T2','Density');
% DOWNSAMPLING a 256-by-256 input: 1 -> 256-by-256, 2 -> 128-by-128, 4 -> 64-by-64.
pixel_stride = 1;       % Keep every pixel_stride-th row/column; no interpolation.
% Optional direct division is shown after loading the maps below.
lower_bounds = [100 5]; % [T1 T2] in ms. This legacy BLIP has no density bounds.
upper_bounds = [5500 550];
% If you divide T1 below, adjust its bounds and dictionary_spacing manually.
dictionary_spacing = [100 10]; % [T1 T2] spacing in ms; smaller spacing costs more time/memory.

%% 1. Settings: change these to try a different experiment
L = 100;                % Number of RF pulses / measured time frames (1 to 100).
r = 16;                % Must divide the reduced row count; r = 1 works for any size.
noise_std = 4;          % Standard deviation of EACH real/imaginary noise part.
% Uniform dictionary grids follow your bounds, including both endpoints.
T1set = unique([lower_bounds(1):dictionary_spacing(1):upper_bounds(1), upper_bounds(1)]);
T2set = unique([lower_bounds(2):dictionary_spacing(2):upper_bounds(2), upper_bounds(2)]);
repetition_time_ms = 10; % Constant illustrative MRF-scale TR; simulator TE = TR/2.
% Original bSSFP-MRF used variable TR between 10.5 and 14 ms:
% Ma et al. (2013), https://doi.org/10.1038/nature11971.
rng(1);                % Fixed seed makes repeated experiments reproducible.
% For a basic noiseless check, set r = 1 and noise_std = 0.
% Dictionary discretization can still cause error in that case.

%% 2. Ground truth and pulse sequence
maps = load(input_file, 'T1', 'T2', 'Density');
% Convert numeric storage classes to double; all maps must describe the same grid.
for name = {'T1','T2','Density'}
    validateattributes(maps.(name{1}), {'numeric'}, {'real','finite','nonempty','2d'});
    maps.(name{1}) = double(maps.(name{1}));
end
assert(isequal(size(maps.T1),size(maps.T2),size(maps.Density)), ...
    'T1, T2 and Density must have the same size.');
assert(all(maps.Density(:)>=0),'Density must be nonnegative.');
% OPTIONAL: uncomment this line to use T1/10 as the new ground truth.
% Adjust the settings above manually. Use the same form for T2 or Density if wanted.
% maps.T1 = maps.T1 / 10;
% Index rows and columns independently, so rectangular images also work.
T1 = maps.T1(1:pixel_stride:end, 1:pixel_stride:end);
T2 = maps.T2(1:pixel_stride:end, 1:pixel_stride:end);
Density = maps.Density(1:pixel_stride:end, 1:pixel_stride:end);
mask = Density > 0;     % Known object support, also passed to the legacy BLIP solver.
[n1, n2] = size(T1);
assert(any(mask(:)),'The reduced input must contain positive Density values.');
% assert(condition,message) stops when the condition is false. Here it
% catches a dictionary that cannot represent the foreground relaxation times.
assert(min(T1(mask)) >= min(T1set) && max(T1(mask)) <= max(T1set) ...
    && min(T2(mask)) >= min(T2set) && max(T2(mask)) <= max(T2set), ...
    'Input times are outside the dictionary range. Edit lower_bounds/upper_bounds.');
assert(r >= 1 && r == floor(r) && mod(n1, r) == 0, ...
    'r must be a positive integer dividing the number of image rows.');

sequence = load(fullfile(project_root, 'data', 'sequence_parameter_middle.mat'), 'RF');
% RF contains flip angles/phases. Only RF is loaded; this demo supplies its own TR.
RF = sequence.RF(1:2:2*L); % Select L pulses, taking every second stored pulse.
TR = repetition_time_ms * ones(1, L); % Shared by data and dictionary simulation.
df = 0;                  % No off-resonance; the current simulator ignores df.

%% 3. Simulate magnetization, Fourier sampling, and measurement noise
% The same Bloch simulator is used here and for the BLIP dictionary below.
bloch = Blochmap_vectorized(RF, TR, T1, T2, df);
% The unit-density Bloch response is n1-by-n2-by-L; density broadcasts over L.
X_clean = Density .* bloch;
zeta = mod(0:L-1, r) + 1;  % Shift the sampled rows by one at each time frame.
noisy_data = complex(zeros(n1, n2, L));
% Zero-filled storage does not mean that unmeasured Fourier entries were acquired.
for t = 1:L
    clean_kspace = fft2(X_clean(:, :, t)) / sqrt(n1*n2);
    % The 1/sqrt(pixel count) factor makes the full Fourier transform unitary.
    sampled_rows = zeta(t):r:n1;
    % Independent Gaussian noise on the acquired Fourier coefficients.
    % Each complex sample has expected squared noise magnitude 2*noise_std^2.
    noise = noise_std * (randn(numel(sampled_rows), n2) ...
                       + 1i*randn(numel(sampled_rows), n2));
    noisy_data(sampled_rows, :, t) = clean_kspace(sampled_rows, :) + noise;
end
% Unmeasured Fourier coefficients remain zero. The spatial mask is not a
% Fourier sampling mask: it is used only as known support in reconstruction.

%% 4. Build the dictionary and reconstruct with the existing BLIP algorithm
% BLIP expects T2 to vary fastest as we move down the dictionary rows.
[T2_grid, T1_grid] = ndgrid(T2set, T1set);
% Flatten the parameter grid for one Bloch simulation of all candidate atoms.
signals = Blochmap_vectorized(RF, TR, T1_grid(:), T2_grid(:), df);
dict = reshape(signals, numel(T1_grid), L);
% Keep raw signal amplitudes: BLIP uses their norms to estimate density.
% BLIP.m currently performs 20 trial iterations, including rejected trials.
fprintf('BLIP: %d-by-%d pixels, %d frames, %d dictionary entries, %dx undersampling.\n', ...
    n1, n2, L, size(dict, 1), r);
tic;
% The legacy helper returns three separate physical maps, not a history structure.
[T1_blip, T2_blip, Density_blip] = BLIP( ...
    L, r, n1, n2, noisy_data, dict, zeta, T1set, T2set, mask);
reconstruction_seconds = toc;

%% 5. Evaluate errors against ground truth
% Relative error (%) = 100 * ||reconstruction - truth||_F / ||truth||_F.
% RMSE is evaluated inside the known object: sqrt(mean(error(mask).^2)).
% Its units are ms for T1/T2 and the original density units for Density.
truth_maps = {T1, T2, Density};
reconstructed_maps = {T1_blip, T2_blip, Density_blip};
map_names = {'T1'; 'T2'; 'Density'};
relative_error_percent = zeros(3, 1);
rmse_in_object = zeros(3, 1);
for j = 1:3
    % This evaluation loop uses truth only after reconstruction has finished.
    truth = truth_maps{j};
    difference = reconstructed_maps{j} - truth;
    relative_error_percent(j) = 100 * norm(difference, 'fro') / norm(truth, 'fro');
    rmse_in_object(j) = sqrt(mean(difference(mask).^2));
end
error_table = table(relative_error_percent, rmse_in_object, ...
    'VariableNames', {'RelativeError_percent', 'RMSE_in_object'}, ...
    'RowNames', map_names);
disp(error_table);
fprintf('BLIP reconstruction time: %.2f seconds.\n', reconstruction_seconds);

%% 6. Plot ground truth, reconstruction, and absolute error for each map
% Ground truth and reconstruction share the ground-truth range in each row.
% Outliers saturate this display scale; reconstruction values are not clipped.
% Absolute-error plots show the full errors in the same physical units.
figure('Name', 'BLIP reconstruction', 'Color', 'w', 'Position', [100 100 1100 800]);
% Rows select T1/T2/density; columns show truth, estimate and absolute error.
tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
plot_names = {'T1 (ms)', 'T2 (ms)', 'Density'};
for j = 1:3
    truth = truth_maps{j};
    reconstruction = reconstructed_maps{j};
    limits = [0, max(truth(:))];
    % Apply the same intensity range to truth and estimate for this parameter.

    ax = nexttile;
    imagesc(truth, limits); axis image off; colorbar; colormap(ax, gray);
    title([plot_names{j}, ': ground truth']);

    ax = nexttile;
    imagesc(reconstruction, limits); axis image off; colorbar; colormap(ax, gray);
    title([plot_names{j}, ': BLIP']);

    ax = nexttile;
    imagesc(abs(reconstruction - truth)); axis image off; colorbar; colormap(ax, hot);
    title(sprintf('Absolute error (relative: %.1f%%)', relative_error_percent(j)));
end
sgtitle(sprintf('BLIP: %d x %d pixels, L = %d, TR = %g ms, r = %d, noise std = %g', ...
    n1, n2, L, repetition_time_ms, r, noise_std));
