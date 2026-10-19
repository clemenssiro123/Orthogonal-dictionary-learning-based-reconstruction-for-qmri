function fig = plot_dictionary_comparison(histories,method_names,options,design,page)
%PLOT_DICTIONARY_COMPARISON Two-method DL progress or body errors.
% FIG = PLOT_DICTIONARY_COMPARISON(HISTORIES,METHOD_NAMES,OPTIONS,DESIGN,PAGE)
% HISTORIES is a two-element cell of nested_dictionary_reconstruction
% histories. METHOD_NAMES and DESIGN colors/styles must use the same order.
% OPTIONS supplies the shared StepNorm, MapTolerance and
% DictionaryTolerance used by both runs. DESIGN is a complete display
% structure including MapNames; partial saved-design overrides are not used.
% PAGE defaults to 'convergence': objective versus steps/time, absolute
% map/D-C changes, and summed inner work. 'errors' uses the three saved
% monitor columns [T1,T2,Density] to plot relative body errors (%).
% State histories include iteration zero; change/work histories start at one.
% Full objectives are comparable only when both methods use the same prior
% weights. No full sequence of intermediate maps is needed for these plots.
% FIG remains editable; no files are written and no solver is run.
% MATLAB R2021b.
% See also PLOT_CONVERGENCE, PLOT_NESTED_CONVERGENCE.
    if nargin<5, page='convergence'; end
    page=validatestring(page,{'convergence','errors'});
    colors=design.MethodColors; labels=design.MapNames;
    switch page
        case 'convergence'
            fig=reconstruction_report_figure('Dictionary comparison: convergence',design);
            tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
            nexttile; hold on;
            for method=1:2
                state=histories{method};
                plot(0:state.iterations,max(state.objective,eps),'Color',colors(method,:),'LineStyle',design.LineStyles{method},'LineWidth',design.LineWidth);
            end
            set(gca,'YScale','log'); grid on; xlabel('Outer step'); title('Full objective');
            legend(method_names,'Location','best');
            nexttile; hold on;
            for method=1:2
                state=histories{method};
                plot(state.time,max(state.objective,eps),'Color',colors(method,:),'LineStyle',design.LineStyles{method},'LineWidth',design.LineWidth);
            end
            set(gca,'YScale','log'); grid on; xlabel('Elapsed seconds'); title('Full objective versus time');
            for change_panel=1:2
                nexttile; hold on;
                for method=1:2
                    state=histories{method}; values=state.step_norm;
                    if change_panel==2, values=state.dictionary_step_norm; end
                    plot(1:state.iterations,max(values,eps),'Color',colors(method,:), ...
                        'LineStyle',design.LineStyles{method},'LineWidth',design.LineWidth);
                end
                tolerance=options.MapTolerance;
                if change_panel==2, tolerance=options.DictionaryTolerance; end
                if tolerance>0, yline(tolerance,':k','Stopping tolerance','HandleVisibility','off'); end
                set(gca,'YScale','log'); grid on; xlabel('Completed outer step'); ylabel('Absolute change');
                if change_panel==1
                    title(sprintf('Successive map change (%s norm)',options.StepNorm));
                else
                    title('Successive D,C change (combined Frobenius norm)');
                end
            end
            nexttile([1 2]); hold on;
            for method=1:2
                state=histories{method};
                stairs(1:state.iterations,sum(state.inner_iterations,2),'Color',colors(method,:),'LineStyle',design.LineStyles{method},'LineWidth',design.LineWidth);
            end
            grid on; xlabel('Completed outer step'); ylabel('D/C update pairs, summed over maps');
            title('Dictionary work per outer step'); % One-step uses 3 pairs: one for each active map.
        case 'errors'
            fig=reconstruction_report_figure('Dictionary comparison: errors',design);
            tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
            for row=1:2
                for j=1:3
                    nexttile; hold on;
                    for method=1:2
                        state=histories{method};
                        horizontal=0:state.iterations;
                        if row==2, horizontal=state.time; end
                        plot(horizontal,state.monitor(:,j),'Color',colors(method,:),'LineStyle',design.LineStyles{method},'LineWidth',design.LineWidth);
                    end
                    grid on; ylabel('Relative error inside body (%)'); title(labels{j});
                    if row==1, xlabel('Outer step'); else, xlabel('Elapsed seconds'); end
                end
            end
            legend(method_names,'Location','best');
    end
end
