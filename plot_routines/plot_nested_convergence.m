function fig = plot_nested_convergence(info,options,design)
%PLOT_NESTED_CONVERGENCE Nested DL objective, changes and inner accuracy.
% FIG = PLOT_NESTED_CONVERGENCE(INFO,OPTIONS,DESIGN)
% INFO is the history returned by nested_dictionary_reconstruction.
% OPTIONS are the options used for that solve, including StepNorm,
% MapTolerance and DictionaryTolerance; they are not inferred from INFO.
% For saved comparisons, pass s.saved_histories{4}, s.saved_options{4},
% and s.design respectively. DESIGN is a complete display structure.
%
% The six panels show objective contributions, absolute map change,
% combined outer D/C change, inner update counts, last inner D/C change
% versus eta_k, and unmasked background density RMS. Last inner change
% differs from the total D/C change across the whole outer update.
% State curves include iteration zero; step curves start at update one.
% FIG remains editable; no files are written and no solver is run.
% MATLAB R2021b.
% See also PLOT_CONVERGENCE, PLOT_DICTIONARY_COMPARISON.
    k=0:info.iterations;
    fig=reconstruction_report_figure('Nested dictionary convergence',design);
    tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
    nexttile;
    colororder(gca,lines(5));
    semilogy(k,max([info.objective info.data info.patch info.spatial info.sparsity],eps),'LineWidth',design.LineWidth);
    grid on; xlabel('Outer step'); title('Objective contributions');
    legend('Total','Data','Patch fit','Spatial','L1','Location','best');
    nexttile;
    semilogy(1:info.iterations,max(info.step_norm,eps),'LineWidth',design.LineWidth);
    if options.MapTolerance>0, yline(options.MapTolerance,':k','Stopping tolerance'); end
    grid on; xlabel('Completed outer step'); ylabel('Absolute change');
    title(sprintf('Successive map change (%s norm)',options.StepNorm));
    nexttile;
    % Outer stopping uses the total D/C change across a complete outer update.
    semilogy(1:info.iterations,max(info.dictionary_step_norm,eps),'LineWidth',design.LineWidth);
    if options.DictionaryTolerance>0, yline(options.DictionaryTolerance,':k','Stopping tolerance'); end
    grid on; xlabel('Completed outer step'); ylabel('Combined Frobenius norm'); title('Successive D,C change');
    nexttile;
    colororder(gca,lines(3));
    plot(1:info.iterations,info.inner_iterations,'LineWidth',design.LineWidth);
    grid on; xlabel('Outer step'); title('Inner D/C update pairs'); legend('T1','T2','Density');
    nexttile;
    colororder(gca,lines(3));
    % Inner accuracy instead tests only the LAST inner change against eta_k.
    semilogy(1:info.iterations,max(info.inner_step_norm,eps),'LineWidth',design.LineWidth); hold on;
    semilogy(1:info.iterations,info.inner_tolerance,'k--','LineWidth',design.LineWidth);
    grid on; xlabel('Outer step'); title('Last inner change versus eta_k');
    legend('T1','T2','Density','eta_k','Location','best');
    nexttile;
    plot(k,info.background_rms(:,3),'LineWidth',design.LineWidth);
    grid on; xlabel('Outer step'); title('Background RMS density');
end
