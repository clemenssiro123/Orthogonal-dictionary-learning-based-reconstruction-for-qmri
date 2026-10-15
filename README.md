# Orthogonal Dictionary Learning for Quantitative MRI

This repository provides a **small-scale MATLAB implementation** of the orthogonal dictionary-learning approach to quantitative MRI (qMRI) presented by [Dong, Hintermüller, and Sirotenko (2025)](https://iopscience.iop.org/article/10.1088/1361-6420/adef74). The emphasis is on small-scale numerical experiments, with adaptations to parameter scaling, initialization, and solver settings.

The examples recover longitudinal relaxation time $T_1$, transverse relaxation time $T_2$, and proton density $\rho$ from simulated, undersampled Fourier measurements generated with a Bloch model.

The main demonstration compares BLIP, projected Levenberg-Marquardt (LM), one-step dictionary learning, and nested dictionary learning using the same measurements and coarse-BLIP initialization.

## Tested environment

The code has been tested on a **Dell XPS 17 laptop running Windows with MATLAB R2021b**.

- MATLAB R2021b.
- Optimization Toolbox for the dictionary-regularized reconstruction routines, which use `quadprog`.

The current implementation is best used at **128 x 128 pixels**, which is the recommended starting resolution for the supplied experiments.

## Quick start

1. Download or clone this repository and open its root folder in MATLAB.
2. Open [`demo_compare_all_reconstructions.m`](demo_compare_all_reconstructions.m).
3. Keep `pixel_stride = 2` to reduce the supplied 256 x 256 input maps to 128 x 128 pixels. Review the input, algorithm, and plotting settings at the top of the script.
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

Adjust `design.RowPanelSize`, `design.RowPanelGap`, and `design.RowShowTitles` in the demo to control the compact layout; set `RowShowTitles=false` to omit headings. `MapLimits`, `ErrorLimits`, `MapColors`, `ErrorColors`, `ZeroBackground`, and `ZoomBox` apply to both the existing comparison grid and the compact PDFs.

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

## Expected runtime

These are **provisional estimates for the entire `demo_compare_all_reconstructions` script**, including all four reconstruction methods, on the Dell XPS 17:

| Reconstruction size | Setting for the supplied input | Estimated wall-clock time |
| --- | --- | --- |
| **128 x 128 (recommended)** | `pixel_stride = 2` | Approximately **20 minutes** |
| 256 x 256 | `pixel_stride = 1` | Approximately **2 hours** |

These values reflect the maintainer's recollection and estimates; they have **not been remeasured for the current version**. Actual runtime depends on hardware, acquisition length, patch size, iteration limits, solver tolerances, and plotting. They should be treated as planning estimates rather than benchmark results.

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

| Script | Purpose |
| --- | --- |
| [`demo_compare_all_reconstructions.m`](demo_compare_all_reconstructions.m) | Compare all four reconstruction methods and export reports. |
| [`demo_compare_dictionary_reconstruction.m`](demo_compare_dictionary_reconstruction.m) | Compare one-step and nested dictionary learning. |
| [`demo_nested_dictionary_reconstruction.m`](demo_nested_dictionary_reconstruction.m) | Run nested dictionary-regularized qMRI reconstruction. |
| [`demo_dictionary_reconstruction.m`](demo_dictionary_reconstruction.m) | Demonstrate orthogonal dictionary learning and denoising on image patches. |
| [`demo_projected_lm_reconstruction.m`](demo_projected_lm_reconstruction.m) | Run the projected LM reconstruction baseline. |
| [`demo_blip.m`](demo_blip.m) | Run the BLIP reconstruction baseline. |

## Input data and outputs

The supplied examples use `data/PWbrain256.mat` for reference parameter maps and `data/sequence_parameter_middle.mat` for the RF pulse sequence. Measurements are simulated from these inputs.

To use your own reference maps, set `input_file` in the selected demo to a MAT file containing equally sized, finite, two-dimensional arrays named `T1`, `T2`, and `Density`. Relaxation times are in milliseconds; density must be nonnegative, with zero marking the background. Review the physical bounds, scaling, and sampling settings for the new input. The demos also include an optional very simple artificial-phantom generator.

The solvers return reconstructed maps and iteration histories. The comparison demo additionally records errors, residuals, timings, and stopping reasons in `results.mat`.

## Reference

Guozhi Dong, Michael Hintermüller, and Clemens Sirotenko. **Dictionary learning based regularization in quantitative MRI: a nested alternating optimization framework.** *Inverse Problems*, **41** (2025), 085007.

[Read the journal article](https://iopscience.iop.org/article/10.1088/1361-6420/adef74) · [DOI: 10.1088/1361-6420/adef74](https://doi.org/10.1088/1361-6420/adef74)

Please cite this article when using the dictionary-learning framework in research.
