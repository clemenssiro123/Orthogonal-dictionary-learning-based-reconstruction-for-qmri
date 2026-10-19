function fig = plot_dictionary_learning_convergence(fixed_info,info,tolerance,design)
%PLOT_DICTIONARY_LEARNING_CONVERGENCE D/C changes across two learning phases.
% FIG = PLOT_DICTIONARY_LEARNING_CONVERGENCE(FIXED_INFO,INFO,TOLERANCE,DESIGN)
% FIXED_INFO and INFO are histories returned by dictionary_reconstruction
% for the fixed-step phase and its warm-started continuation respectively.
% Their step_norm vectors are joined in execution order; each value is a
% combined dictionary/coefficient iterate change, starting at update one.
% TOLERANCE draws the horizontal stopping line. The vertical line at
% FIXED_INFO.iterations marks where continuation begins.
% DESIGN is a complete display structure with at least three method colors,
% as used by demo_dictionary_reconstruction; color three draws the curve.
% FIG remains editable; no files are written and no solver is run.
% MATLAB R2021b.

    fig=reconstruction_report_figure('Dictionary convergence',design);
    semilogy([fixed_info.step_norm; info.step_norm],'Color',design.MethodColors(3,:), ...
        'LineWidth',design.LineWidth);
    % Join the histories; distinguish the continuation point from the tolerance.
    hold on;
    yline(tolerance,'--','Tolerance');
    xline(fixed_info.iterations,':','Continuation');
    xlabel('D/C update pairs'); ylabel('Combined iterate change'); grid on;
    title('Paper stopping criterion');
end
