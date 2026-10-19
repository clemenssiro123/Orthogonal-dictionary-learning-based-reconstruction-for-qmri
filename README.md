# Orthogonal Dictionary Learning for Quantitative MRI

This repository provides a **small-scale MATLAB implementation** of the orthogonal dictionary-learning approach to quantitative MRI (qMRI) presented by [Dong, Hintermüller, and Sirotenko (2025)](https://iopscience.iop.org/article/10.1088/1361-6420/adef74). The emphasis is on small-scale numerical experiments, with adaptations to parameter scaling, initialization, and solver settings.

The examples recover longitudinal relaxation time $T_1$, transverse relaxation time $T_2$, and proton density $\rho$ from simulated, undersampled Fourier measurements generated with a Bloch model.

The main demonstration compares BLIP, projected Levenberg-Marquardt (LM), one-step dictionary learning, and nested dictionary learning using the same measurements and coarse-BLIP initialization.

## Tested environment

The code has been tested on a **Dell XPS 17 laptop running Windows with MATLAB R2021b**.
The notebook is an XPS 17 9710 with an Intel Core i7-11800H processor (8 cores, 16 threads), 64 GB RAM, and an NVIDIA GeForce RTX 3060 Laptop GPU, running Windows 11 Home.

- MATLAB R2021b.
- Optimization Toolbox for the dictionary-regularized reconstruction routines, which use `quadprog`.

The current implementation is best used at **128 x 128 pixels**, which is the recommended starting resolution for the supplied experiments.

## Quick start

1. Download or clone this repository and open its root folder in MATLAB.
2. Open [`demo_compare_all_reconstructions.m`](demo_compare_all_reconstructions.m).
3. The script loads the saved **128 x 128 baseline** by default (`pixel_stride = 2` on a 256 x 256 source). Select `reconstruction_parameters_256()` in its parameter section for full resolution. Review the input and plotting settings in the demo and the algorithm settings in the selected preset file.
4. Run the script from the MATLAB Command Window:

```matlab
demo_compare_all_reconstructions
```

The script adds its helper folders to the MATLAB path automatically. It displays reconstructions and diagnostics, and saves `results.mat`, four PDF reports (`maps.pdf`, `errors.pdf`, `convergence.pdf`, and `profiles_and_summary.pdf`), and eight compact PDFs in a new `data/comparison_all_<timestamp>/` folder. Generated results are excluded from Git.

The eight additional PDFs have small margins and linked zoom insets:

- `T1_map_comparison.pdf` and `T1_error_comparison.pdf`
- `T2_map_comparison.pdf` and `T2_error_comparison.pdf`
- `Density_map_comparison.pdf` and `Density_error_comparison.pdf`
- `initialization_maps.pdf` and `initialization_errors.pdf`

Each of the six comparison files contains one row in this order: **ground truth, BLIP, Levenberg-Marquardt, one-step dictionary learning, nested dictionary learning**. Error panels show absolute differences from ground truth (milliseconds for T1/T2, original density units); the ground-truth error panel is zero. Every row shares one colorbar. Each initialization PDF contains **T1, T2, and density**, with a separate colorbar for each parameter. Color limits and the zoom region match the comparison PDFs.

Adjust `design.RowPanelSize`, `design.RowPanelGap`, and `design.RowShowTitles` in the demo to control the compact layout; set `RowShowTitles=false` to omit comparison method headings. Both initialization PDFs always omit headings and parameter labels, keeping the colorbars and numeric ticks. `MapLimits`, `ErrorLimits`, `MapColors`, `ErrorColors`, `ZeroBackground`, and `ZoomBox` apply to both the existing comparison grid and the compact PDFs.

The eight compact PDFs are exported without opening figure windows by default (`design.ShowRowFigures=false`); each hidden figure is closed after export. Set `design.ShowRowFigures=true` to inspect those figures on screen as well. The master switch `design.ShowFigures=false` hides all report and compact figures while still saving the PDFs.

To generate only the eight compact PDFs from an existing experiment, without rerunning the solvers:

```matlab
run_folder = fullfile(pwd,'data','comparison_all_<timestamp>'); % Choose your saved run.
s = load(fullfile(run_folder,'results.mat'));
s.design.ShowFigures = false;
[pdf_paths, zoom_box] = export_reconstruction_rows(s.truth,s.reconstructions, ...
    s.u0,s.mask,s.method_names,s.design,run_folder);
```

This replaces any compact PDFs with the same filenames; the saved `results.mat` and four existing reports are unchanged. Missing compact-layout settings in older runs use the defaults.

For the supplied input, `pixel_stride = 1` uses the full 256 x 256 resolution; this requires substantially more computation. The undersampling factor `r` must divide the reconstructed image's row count.

## Saved parameter sets

Two editable MATLAB functions store the acquisition settings, physical bounds/scales, BLIP grids and initialization options, and the BLIP, LM, one-step DL and nested DL solver options:

- [`reconstruction_parameters_128.m`](reconstruction_parameters_128.m): the current four-method comparison settings, saved as the 128 x 128 baseline.
- [`reconstruction_parameters_256.m`](reconstruction_parameters_256.m): the standard 256 x 256 configuration from the supplied four-method comparison, with `pixel_stride = 1`, `r = 32`, noise SD `1`, and 80 outer steps for both DL methods. It is independent of the 128 x 128 preset.

Load either set as a MATLAB structure:

```matlab
parameters = reconstruction_parameters_128();
% parameters = reconstruction_parameters_256();
```

In `demo_compare_all_reconstructions.m`, select the desired function in **Load saved experiment and algorithm parameters**. These functions return a plain MATLAB structure without running a reconstruction or changing files. Keeping configuration in small `.m` functions makes the values readable, records edits in Git, and allows derived quantities such as dictionary grids and initial damping to be computed from the same settings.

### Changing parameters

**To save a new default, edit the values inside the selected preset function.** For example, to tune the 256 x 256 configuration:

1. Open `reconstruction_parameters_256.m` in the MATLAB Editor (`edit reconstruction_parameters_256`).
2. Change acquisition variables such as `r` or `noise_std`, or values in the relevant options structure. For example, change `'MaxIterations',15` to `'MaxIterations',20` inside `lm_options=struct(...)` to increase the LM budget. This is an editing example, not a recommendation to use 20 iterations.
3. Save the file. The final `parameters=struct();` block only collects the settings for return; normally leave that block unchanged.
4. Select `parameters=reconstruction_parameters_256();` in the demo and run it. The 128 x 128 preset remains unchanged.

Within each preset, `nested_options=one_options` copies the one-step DL settings before adding the nested-specific controls. Editing the original `one_options` definition therefore changes both methods' shared settings. To change only nested DL, edit or add a `nested_options.<field>=...` assignment after that copy and before the final return block. To change only one-step DL, place its override after the copy as well.

**For a temporary experiment, override the returned structure in the demo.** Place overrides immediately after the preset-loading line and before the demo copies its fields into working variables:

```matlab
parameters = reconstruction_parameters_256();
parameters.lm_options.MaxIterations = 20;       % Example one-run override.
parameters.nested_options.InnerMaxIterations = 1500;
```

These assignments leave the preset file unchanged and apply on every run until removed from the demo. After loading, `parameters.one_options` and `parameters.nested_options` are independent copies; change both explicitly if an override should affect both. Change dependent settings such as bounds, grid spacing or `r` inside the preset function so the grids, BLIP step sizes and LM damping are rebuilt together. Assigning `parameters.r` after loading does not recompute those derived values.

### Current choices and their roles

The 128 x 128 preset records the current baseline, while the 256 x 256 preset records the supplied standard configuration. Their weights, iteration budgets and tolerances are practical choices for these setups; they are not claimed to be optimal for other inputs, noise levels or resolutions. Values below apply to both presets unless a resolution is specified.

| Setting in the preset | Current value | Purpose and effect of changing it |
| --- | --- | --- |
| `L`, `r`, `TR_ms`, `noise_std`, `seed` | 128: `100`, `16`, `12`, `1.5`, `1`; 256: `100`, `32`, `12`, `1`, `1` | Define the common simulated acquisition. Noise SD applies separately to real and imaginary parts. The seed fixes the random draw for a given setup. |
| `lower_bounds`, `upper_bounds`, `scales`, `reference_times` | `[100 10 0]`, `[5500 600 110]`, `[5000 500 100]`, `[1000 100]` | Define the admissible physical maps and shifted, scaled optimization coordinates. Time values are in ms. Changing scales also changes the DL prior in physical units, so review its weights. |
| `coarse_spacing`, `initial_options`, `initial_use_body_mask` | `[200 50]` ms; 29 x 14 grid; at most 20 trials; mask enabled | Build a common coarse-BLIP start at lower dictionary cost. Rejected backtracking trials count toward this initialization budget. The mask zeros initial background density only. |
| `fine_spacing`, `blip_options.MaxIterations` | `[100 10]` ms; 55 x 60 grid; at most 30 accepted updates | Refine the BLIP time grid and allow further reconstruction updates. Smaller spacing increases dictionary size and matching cost. The coarse grid is included in the fine dictionary. |
| `lm_options` damping and budget | `InitialDamping=r^2` (256: `1024`), `DampingDecay=0.15`, `MinDamping=1e-8`, `MaxIterations=15` | Control stabilization and the number of projected LM updates. A smaller decay reduces damping faster; more updates can improve the data fit but also fit more noise. |
| DL `Alpha`, `Lambda`, `Beta`, `PatchSize` | Both: `Alpha=0.7`, `[8 8]` patches; 128: `Lambda=[10.75 37.5 37.5]*1.5`, `Beta=[.009 .009 .009]*1.5`; 256: `Lambda=[10.125 56.25 56.25]`, `Beta=[.0126 .0126 .0126]` | Weight spatial smoothness, patch fit and coefficient sparsity on the scaled maps. The coefficient penalty uses `Lambda.*Beta`; all three-component vectors follow `[T1 T2 Density]`. |
| DL `MaxIterations`, `InnerStopRule` | 128: 100 outer updates; 256: 80; `'one-step'` or `'tolerance'` | Set the outer budget and amount of dictionary learning between map updates. Each preset shares the prior weights across both DL methods. |
| `h`, `GaussNewton`, DL `StepNorm` and QP tolerances | `1`; `'approximate'` for LM and DL; `'H01'`; optimality/function/PCG: `[1e-7 0 1e-6]` | Set grid spacing, curvature approximation, DL step norm and map-subproblem accuracy. |
| `nested_options.InnerMaxIterations`, `InnerToleranceScale`, `InnerDecay` | `1000`, `0.005`, `0.75` | Cap inner dictionary updates and set their decreasing change tolerance. A smaller scale requests a tighter inner solve, potentially increasing cost or reaching the cap. |

The LM choice `InitialDamping=r^2` follows the numerical choice in [Dong et al., Section 4.4](https://arxiv.org/pdf/1903.11979), where their undersampling factor `s` corresponds to this code's `r`. The decay `0.15` and 15-update budget are settings for this implementation, not values prescribed by that paper. LM has no spatial or dictionary penalty in this comparison.

`NoiseLevel=[]` disables discrepancy stopping for BLIP and LM. To try it in the four-method demo, uncomment the `lm_options.NoiseLevel=expected_noise_norm` and/or `blip_options.NoiseLevel=expected_noise_norm` lines after data generation. `NoiseLevel` is a noise **norm**, not `noise_std`; `DiscrepancyFactor` sets the residual threshold multiplier. Iteration limits are caps, and the saved `stop_reason` records the actual termination. LM's `Tolerance` tests its projected-gradient measure, while DL's `MapTolerance` and `DictionaryTolerance` test outer changes; `QPOptimalityTolerance` concerns a map subproblem. Meeting a step tolerance alone does not quantify closeness to stationarity without the constants in the corresponding bounds.

The DL outer change test follows the two-block structure of [Algorithms 2 and 4 in the paper](https://arxiv.org/pdf/2506.11977): both the absolute change between successive scaled maps in the selected `StepNorm` and the combined Frobenius change of `(D,C)` must meet their tolerances. The saved presets use `H01`; `H1` also includes the map-value term. The local prox-linear solver uses `StepTolerance` for its absolute map change in `StepNorm` (default `1e-6`), with `Tolerance` accepted as a legacy alias when `StepTolerance` is omitted; zero tests exact zero change. Its `step_norm_absolute` history records the stopping quantity, while `step_norm` retains the relative change as an auxiliary diagnostic. Prox-linear and DL do not compute a projected-gradient measure. LM retains its projected-gradient stopping rule and history.

### Loading presets in other demos

The standalone LM, nested DL and two-method DL comparison demos have an **Optional saved qMRI parameters** section after their manual settings. Leave `parameters = []` to keep their existing defaults, or uncomment one preset-loading line. The two-method comparison uses `nested_options` for both methods and changes only `InnerStopRule`; use the four-method comparison to test independently configured one-step and nested options. The legacy BLIP demo can load the acquisition and dictionary grid, but retains its own fixed 20-trial budget and support-mask behavior. The image-denoising demo shows how to reuse only `pixel_stride`, since qMRI solver parameters do not apply to that problem.

The preset sizes assume a **256 x 256 source image**. Input-file selection, artificial-phantom generation and plotting stay in each demo; standalone demos also retain their own initialization. Loading a preset in the Command Window alone does not configure a demo because demos start with `clearvars`: select it inside the script. The four-method demo continues to save the actual experiment settings and solver options in `results.mat`.

## Expected runtime

The 128 x 128 estimate is provisional for the entire `demo_compare_all_reconstructions` script on the Dell XPS 17. The 256 x 256 timing comes from the supplied single-run comparison using the standard preset:

| Reconstruction size | Setting for the supplied input | Runtime information |
| --- | --- | --- |
| **128 x 128 (recommended)** | `pixel_stride = 2` | Approximately **20 minutes** |
| 256 x 256 | `pixel_stride = 1` | Approximately **5 hours 33 minutes** for all four solves plus one common initialization (19,963.10 s) |

The 128 x 128 estimate reflects the maintainer's recollection and has not been remeasured for the current version. The reported 256 x 256 solver times include monitoring/setup; simulation and report generation add to the full script runtime. All four methods in that run stopped at `max_iterations`, which does not certify convergence. Actual runtime depends on hardware, acquisition length, patch size, iteration limits, solver tolerances, and plotting.

The saved comparison table reports `SolveSeconds` for each method and `TotalSeconds` for that method plus the shared initialization. These per-method values do not represent the full script's wall-clock time, which also includes simulation and report generation.

## Mathematical description

For context, the underlying approach is summarized by the paper's discrete problem $(P_2)$. Written with consistent component indices and explicit constraints, it is

```math
\begin{aligned}
\min_{u,D,C}\quad
& \frac{h^2}{2}\bigl\lVert F_d(u)-f^\delta \bigr\rVert_2^2
+ \frac{\alpha}{2}\bigl\lVert \nabla u \bigr\rVert_{U_1}^2
\\[0.6em]
& {}+ \sum_{j=1}^{3}\lambda_j\left(
\frac{1}{2}\bigl\lVert P[u_j/M_j]-D_jC_j \bigr\rVert_F^2
+ \beta_j\bigl\lVert C_j \bigr\rVert_1\right)
\\[0.8em]
\text{subject to}\quad
& u\in U_{\mathrm{ad}},\qquad D_j^{\mathsf T}D_j=I,\qquad j=1,2,3.
\end{aligned}
```

Here $F_d$ combines Bloch dynamics and sampled Fourier measurements, $f^\delta$ denotes noisy data, $P$ extracts image patches, and $M_j$ scales each parameter map. Orthogonal dictionaries $D_j$ and sparse coefficients $C_j$ describe local image structure; $U_{\mathrm{ad}}$ imposes parameter bounds. The gradient term provides spatial regularization.

**Implementation convention:** MATLAB arrays use the channel order `[T1, T2, Density]`. Both patch and spatial penalties act on shifted, scaled maps $x_j=(u_j-u_{\mathrm{ref},j})/s_j$, controlled by `ReferenceTimes` and `Scales`. The algorithm alternates dictionary/coefficient updates with bounded LM/Gauss-Newton map updates. The current demo settings are starting configurations, rather than an exact reproduction of the paper's experiments.

## Available examples

Run any of these scripts from the MATLAB Editor or enter its name without `.m` in the Command Window. Each script clears the workspace and closes existing figures. Sizes below refer to the supplied input and current `pixel_stride` settings.

| Script | Purpose and initialization | Default image size |
| --- | --- | --- |
| [`demo_compare_all_reconstructions.m`](demo_compare_all_reconstructions.m) | Compare BLIP, projected LM, one-step DL and nested DL using one simulated data set and one common coarse-BLIP start. Saves numerical results and comparison reports. | 128 x 128 |
| [`demo_compare_dictionary_reconstruction.m`](demo_compare_dictionary_reconstruction.m) | Compare one-step and nested DL using the same measurements, uniform foreground start, and initial dictionaries/coefficients. Only the inner dictionary stopping rule differs. | 128 x 128 |
| [`demo_nested_dictionary_reconstruction.m`](demo_nested_dictionary_reconstruction.m) | Run dictionary-regularized qMRI reconstruction from a uniform foreground start. Refines dictionaries/coefficients to an inner tolerance before each map update. | 128 x 128 |
| [`demo_dictionary_reconstruction.m`](demo_dictionary_reconstruction.m) | Denoise a single image (`T1` by default) using overlapping patches. Runs a fixed number of dictionary updates, then continues from that result to a tolerance. Does not simulate Fourier measurements or solve the qMRI inverse problem. | 256 x 256 |
| [`demo_projected_lm_reconstruction.m`](demo_projected_lm_reconstruction.m) | Run projected LM from a uniform foreground start, with physical L2 damping and no spatial or dictionary penalty. Plots residuals, successive map changes, damping and projection diagnostics; retains projected-gradient stopping. | 256 x 256 |
| [`demo_blip.m`](demo_blip.m) | Run the legacy `helper_functions/BLIP.m` from zero signal, using a Bloch dictionary and a known support mask for the returned density. Its budget is 20 trials, including rejected trials. | 256 x 256 |

For comparisons between reconstruction methods, use `demo_compare_all_reconstructions`: the standalone demos have different acquisition, noise and initialization settings, so their default results are not directly comparable. The four-method demo uses `blip_reconstruction.m`, which provides explicit bounds, stopping options and iteration histories; the standalone BLIP demo uses the legacy interface. One-step DL is available in both comparison scripts.

The LM defaults also differ: the four-method comparison uses initial damping `r^2`, decay `0.15` and at most 15 iterations; the standalone LM demo uses initial damping `1`, decay `0.6` and at most 14 iterations for its uniform start. Both currently disable discrepancy stopping. These are experiment settings, not universal parameter choices; edit the selected preset's `lm_options`, or the standalone script's manual settings when preset loading is disabled.

The two qMRI dictionary-learning demos and the four-method comparison require Optimization Toolbox. The standalone image-denoising, LM and BLIP demos do not require that toolbox.

## Input data and outputs

The qMRI reconstruction demos use `data/sequence_parameter_middle.mat` for the RF pulse sequence. The four-method comparison currently generates artificial 256 x 256 reference maps and subsamples them according to the selected preset; its input section also offers `data/PWbrain256.mat`, which is the standalone demos' default. Measurements are simulated from the selected maps. The image-denoising demo loads only the selected image variable and adds noise directly to the image; it does not use the RF sequence.

To use your own reference maps in a qMRI demo, set `input_file` to a MAT file containing equally sized, finite, two-dimensional arrays named `T1`, `T2`, and `Density`. Relaxation times are in milliseconds; density must be nonnegative, with zero marking the background. Review the physical bounds, scaling, and sampling settings for the new input. For image denoising, supply just one real, finite, two-dimensional numeric array and set `image_variable` to its name; adjust `image_scale` and `noise_std` to its units. The demos also include an optional very simple artificial-phantom generator.

The demos leave reconstructed maps and available iteration histories in the MATLAB workspace. Only the four-method comparison automatically saves `results.mat`, including errors, residuals, timings, stopping reasons and experiment settings. Save workspace variables explicitly if numerical results from the other demos are needed later.

The other demos use the same plot style: shared map/error color limits, linked zoom insets, compact PDF rows, and profile reports. They export to a new time-stamped folder under `data/demo_reports/`. Edit each demo's `design` settings to choose the zoom, profile row, background display, or figure visibility; the common defaults are in `reconstruction_plot_design.m`. Convergence PDFs show successive-iterate changes in place of projected-gradient stationarity plots. DL reports show separate absolute map and combined `(D,C)` changes, with dotted stopping-tolerance lines. Prox-linear uses its absolute step metric; LM shows the physical L2 change as a diagnostic while retaining projected-gradient stopping. Change histories begin at the first completed update. The denoising demo exports a single image channel, while legacy BLIP omits initialization and convergence reports because its interface does not return those quantities.

PDF export normally uses vector graphics. If MATLAB's vector renderer fails, `export_reconstruction_pdf.m` retries on a fresh hidden figure, then falls back to an image page at `design.Resolution` DPI with a warning if necessary; existing PDF pages are preserved when an export fails.

## Reference

Guozhi Dong, Michael Hintermüller, and Clemens Sirotenko. **Dictionary learning based regularization in quantitative MRI: a nested alternating optimization framework.** *Inverse Problems*, **41** (2025), 085007.

[Read the journal article](https://iopscience.iop.org/article/10.1088/1361-6420/adef74) · [DOI: 10.1088/1361-6420/adef74](https://doi.org/10.1088/1361-6420/adef74)

Please cite this article when using the dictionary-learning framework in research.
