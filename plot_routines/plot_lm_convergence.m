function fig = plot_lm_convergence(info,design,page)
%PLOT_LM_CONVERGENCE LM convergence or solver diagnostics.
% FIG = PLOT_LM_CONVERGENCE(INFO,DESIGN,PAGE)
% INFO is the history returned by projected_lm_reconstruction, including
% INFO.options; use s.saved_histories{2} for a saved four-method run.
% DESIGN is a complete display structure (for example s.design), not a
% partial override. ShowFigures controls visibility.
% PAGE defaults to 'convergence': residual/noise threshold, absolute
% physical-L2 map change, and background density RMS. 'diagnostics' shows
% damping, linear-solve residual, clipped-entry fraction and relative step.
% State histories include update zero; step histories start at update one.
% The solver's projected-gradient stopping rule is unaffected by plotting.
% FIG remains editable; no files are written and no solver is run.
% MATLAB R2021b.
% Example: fig=plot_lm_convergence(s.saved_histories{2},s.design,'diagnostics');
    if nargin<3, page='convergence'; end
    page=validatestring(page,{'convergence','diagnostics'});
    use_discrepancy=~isempty(info.options.NoiseLevel);
    linear_tolerance=info.options.LinearTolerance;
    switch page
        case 'convergence'
            fig=reconstruction_report_figure('Projected LM: convergence',design);
            tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
            nexttile; semilogy(0:info.iterations,max(info.data_residual_norm,realmin),'LineWidth',design.LineWidth);
            if use_discrepancy && info.discrepancy_threshold>0
                hold on; yline(info.discrepancy_threshold,'k--','tau*delta');
            end
            grid on; xlabel('Update'); ylabel('Residual norm'); title('Measured data and noise threshold');
            nexttile; semilogy(1:info.iterations,max(info.step_norm_absolute,realmin),'LineWidth',design.LineWidth);
            grid on; xlabel('Completed update'); ylabel('Absolute physical L2 norm'); title('Successive map change');
            nexttile; plot(0:info.iterations,info.background_rms(:,3),'LineWidth',design.LineWidth);
            grid on; xlabel('Update'); ylabel('RMS density'); title('Density leaking into background');
        case 'diagnostics'
            % Linear accuracy is measured before clipping to the box bounds;
            % the relative step is measured after that projection.
            fig=reconstruction_report_figure('Projected LM: solver diagnostics',design);
            tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
            nexttile; semilogy(1:info.iterations,info.damping,'LineWidth',design.LineWidth);
            grid on; xlabel('Update'); ylabel('Lambda'); title('Scheduled damping');
            nexttile; semilogy(1:info.iterations,max(info.linear_relative_residual,realmin),'LineWidth',design.LineWidth);
            hold on; yline(linear_tolerance,'k--'); grid on;
            xlabel('Update'); ylabel('Relative residual'); title('Linear system before projection');
            nexttile; plot(1:info.iterations,100*info.projection_fraction,'LineWidth',design.LineWidth);
            grid on; xlabel('Update'); ylabel('Percent of entries'); title('Entries clipped to bounds');
            nexttile; semilogy(1:info.iterations,max(info.step_norm,realmin),'LineWidth',design.LineWidth);
            grid on; xlabel('Update'); ylabel('Relative step'); title('Actual step after projection');
    end
end
