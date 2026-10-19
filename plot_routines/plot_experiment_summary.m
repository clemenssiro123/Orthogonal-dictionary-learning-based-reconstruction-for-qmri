function fig = plot_experiment_summary(results,design)
%PLOT_EXPERIMENT_SUMMARY Four-method numerical/settings summary.
% FIG = PLOT_EXPERIMENT_SUMMARY(RESULTS,DESIGN)
% RESULTS is the four-method results structure, MAT filename, or run folder.
% Reads comparison, experiment, initial_history, initial_options and
% saved_options to reproduce the original summary page. In particular,
% tolerances and timings come from this run, not the current preset files.
% SolveSeconds includes each solver's setup/monitoring; TotalSeconds adds
% the common initialization. Errors refer to body pixels; background RMS
% uses the unmasked reconstructions. Stop reasons are reported verbatim.
% Optional DESIGN overrides saved display fields, including report fonts.
% FIG remains editable; no files are written and no solver is run.
% MATLAB R2021b.
% See also EXPORT_COMPARISON_REPORTS.
    if nargin<2, design=[]; end
    [results,design]=comparison_plot_inputs(results,design);
    [n,m,~]=size(results.truth); e=results.experiment; t=results.comparison;
    L=e.Frames; r=e.Undersampling; noise_std=e.NoiseStd; seed=e.Seed; TR_ms=e.TR_ms;
    scales=e.Scales; reference_times=e.ReferenceTimes; h=e.H;
    initial_seconds=e.InitialSeconds; initial_use_body_mask=e.InitialUseBodyMask;
    initial_history=results.initial_history;
    coarse_T1=results.initial_options.T1Grid; coarse_T2=results.initial_options.T2Grid;
    method_names=results.method_names; method_options=results.saved_options;
    blip_options=method_options{1}; lm_options=method_options{2}; nested_options=method_options{4};
    iterations=t.Steps; seconds=t.SolveSeconds; total_seconds=t.TotalSeconds;
    relative_error=t{:,{'T1Error_percent','T2Error_percent','DensityError_percent'}};
    residual_norm=t.DataResidual; background_density=t.BackgroundDensityRMS; stop_reason=t.StopReason;
    fig=comparison_report_figure('Experiment summary',design); ax=axes(fig,'Position',[.04 .04 .92 .9]); axis(ax,'off');
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

end
