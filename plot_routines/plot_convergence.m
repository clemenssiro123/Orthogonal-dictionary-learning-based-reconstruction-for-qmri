function fig = plot_convergence(results,quantity,design)
%PLOT_CONVERGENCE One convergence plot or a full comparison-report page.
% FIG = PLOT_CONVERGENCE(RESULTS,QUANTITY,DESIGN)
% RESULTS is the four-method results structure, MAT filename, or run folder;
% curves use saved_histories, saved_options, comparison and experiment.
% QUANTITY defaults to 'residual'. Single-axes choices:
%   residual, residual-time, objective, map-change, dictionary-change,
%   inner-iterations, qp-optimality, t1-error, t2-error, density-error.
% Residual/error curves compare all four methods. Objective, changes, inner
% work and QP curves compare the two DL methods (entries 3 and 4).
% Full report pages: 'fit' (residuals/timings/background), 'errors' (body
% errors versus steps/time), and 'dictionary' (five DL diagnostics).
%
% State curves include iteration zero. Change/work curves start at update
% one. Map changes are absolute in the configured solver norm; D/C changes
% use the combined Frobenius norm. Dotted lines show saved stopping
% tolerances. QP optimality describes the map subproblem, not outer stopping.
% Objectives are comparable between DL methods only for equal prior weights.
% Error curves use relative body errors (%); residuals use the common
% measured-data norm. Log-plot floors change display values only.
%
% Optional DESIGN overrides saved display fields. FIG remains editable;
% no files are written and no solver is run. MATLAB R2021b.
% Example: fig=plot_convergence(s,'map-change',struct('ShowFigures',true));
% See also PLOT_COMPARISON, EXPORT_COMPARISON_REPORTS.
    %% 1. Read the saved run and choose one plot or a complete report page
    if nargin<2 || isempty(quantity), quantity='residual'; end
    if nargin<3, design=[]; end
    [results,design]=comparison_plot_inputs(results,design);
    quantity=validatestring(quantity,{'residual','residual-time','objective', ...
        'map-change','dictionary-change','inner-iterations','qp-optimality', ...
        't1-error','t2-error','density-error','fit','errors','dictionary'});
    histories=results.saved_histories; method_options=results.saved_options;
    method_names=results.method_names; map_names=design.MapNames;
    noise_norm=results.experiment.NoiseNorm; seconds=results.comparison.SolveSeconds;
    initial_seconds=results.experiment.InitialSeconds;
    background_density=results.comparison.BackgroundDensityRMS;
    %% 2. Recreate report layouts, or delegate to the single-axes plot below
    switch quantity
        case 'fit'
            fig=comparison_report_figure('Data fit and runtime',design); tiledlayout(2,2,'TileSpacing','compact');
            for panel=1:2
                % Use iteration count in panel 1 and elapsed solver time in panel 2.
                nexttile; hold on;
                for method=1:4
                    state=histories{method}; horizontal=0:state.iterations;
                    if panel==2, horizontal=state.time; end
                    draw_reconstruction_curve(horizontal,state.comparison_residual,method,design);
                end
                if noise_norm>0, yline(noise_norm,':k','Actual noise norm','HandleVisibility','off'); end
                set(gca,'YScale','log'); grid on; ylabel('Measured residual norm');
                if panel==1, xlabel('Completed step'); else, xlabel('Solver seconds (including setup)'); end
                title('Same data-fit quantity for all methods'); legend(method_names,'Location','best');
            end
            nexttile; bar([seconds initial_seconds*ones(4,1)],'stacked');
            % The stacked bars distinguish work unique to a method from the common start.
            set(gca,'XTick',1:4,'XTickLabel',method_names); ylabel('Seconds'); grid on;
            legend({'Solver + its setup','Shared coarse initialization'},'Location','best'); title('Single-run time costs');
            nexttile; bar(background_density); set(gca,'XTick',1:4,'XTickLabel',method_names);
            ylabel('RMS density outside body'); grid on; title('Background leakage (unmasked reconstruction)');
            sgtitle('Data fit and computational cost; iteration costs differ across methods');

        case 'errors'
            fig=comparison_report_figure('Errors over iterations and time',design); tiledlayout(2,3,'TileSpacing','compact');
            % MonitorFcn stored three body errors at every state; plot them without rerunning solvers.
            for row=1:2
                for j=1:3
                    nexttile; hold on;
                    for method=1:4
                        state=histories{method}; horizontal=0:state.iterations;
                        if row==2, horizontal=state.time; end
                        draw_reconstruction_curve(horizontal,state.monitor(:,j),method,design);
                    end
                    grid on; title(map_names{j}); ylabel('Relative body error (%)');
                    if row==1, xlabel('Completed step'); else, xlabel('Solver seconds (including setup)'); end
                    if j==1, legend(method_names,'Location','best'); end
                end
            end
            sgtitle('Reconstruction quality, including the common initial maps at step zero');

        case 'dictionary'
            fig=comparison_report_figure('Dictionary convergence',design); tiledlayout(3,2,'TileSpacing','compact');
            % Full objectives are comparable between the DL methods only for equal prior weights.
            % Iterate-change histories start at one; dotted lines mark stopping tolerances.
            for panel=1:5
                if panel==1, nexttile([1 2]); else, nexttile; end
                hold on;
                for method=3:4
                    state=histories{method}; horizontal=0:state.iterations;
                    switch panel
                        case 1, vertical=state.objective; label='Full regularized objective';
                        case 2
                            horizontal=1:state.iterations; vertical=state.step_norm;
                            label='Successive map change (solver norm)';
                        case 3
                            horizontal=1:state.iterations; vertical=state.dictionary_step_norm;
                            label='Successive D,C change (combined Frobenius norm)';
                        case 4
                            horizontal=1:state.iterations; vertical=sum(state.inner_iterations,2);
                            label='D/C pairs per outer step (sum over maps)';
                        case 5
                            horizontal=1:state.iterations; vertical=state.qp_firstorderopt;
                            label='QP first-order optimality';
                    end
                    draw_reconstruction_curve(horizontal,max(vertical,eps),method,design);
                    if panel==2 || panel==3
                        tolerance=method_options{method}.MapTolerance;
                        if panel==3, tolerance=method_options{method}.DictionaryTolerance; end
                        if tolerance>0
                            yline(tolerance,':','Color',design.MethodColors(method,:),'HandleVisibility','off');
                        end
                    end
                end
                set(gca,'YScale','log'); grid on; xlabel('Completed step'); title(label);
                legend(method_names(3:4),'Location','best');
            end
            sgtitle('DL diagnostics: absolute iterate changes; dotted lines mark stopping tolerances');

        otherwise
            fig=single_curve(histories,method_options,method_names,design,quantity,noise_norm);
    end
end

function fig=single_curve(histories,options,names,design,quantity,noise_norm)
% Map the requested quantity to its history field and appropriate methods.
% Four-method saved histories are ordered BLIP, LM, one-step DL, nested DL.
    methods=1:numel(histories); is_step=false; is_log=true; tolerance_field='';
    switch quantity
        case {'residual','residual-time'}
            field='comparison_residual'; label='Measured residual norm';
        case 'objective'
            methods=3:4; field='objective'; label='Full regularized objective';
        case 'map-change'
            methods=3:4; field='step_norm'; is_step=true;
            label='Successive map change (solver norm)'; tolerance_field='MapTolerance';
        case 'dictionary-change'
            methods=3:4; field='dictionary_step_norm'; is_step=true;
            label='Successive D,C change (combined Frobenius norm)'; tolerance_field='DictionaryTolerance';
        case 'inner-iterations'
            methods=3:4; field='inner_iterations'; is_step=true; is_log=false;
            label='D/C pairs per outer step (sum over maps)';
        case 'qp-optimality'
            methods=3:4; field='qp_firstorderopt'; is_step=true; label='QP first-order optimality';
        otherwise
            field='monitor'; is_log=false;
            channel=find(strcmp(quantity,{'t1-error','t2-error','density-error'}));
            label=[design.MapNames{channel} ' relative body error (%)'];
    end
    fig=comparison_report_figure(label,design); hold on;
    for method=methods
        state=histories{method}; vertical=state.(field);
        % States have an initial sample; update-only quantities do not.
        horizontal=0:state.iterations;
        if is_step, horizontal=1:state.iterations; end
        if strcmp(quantity,'residual-time'), horizontal=state.time; end
        if strcmp(field,'inner_iterations'), vertical=sum(vertical,2); end
        if strcmp(field,'monitor'), vertical=vertical(:,channel); end
        % Keep log axes finite without modifying the original saved history.
        if is_log, vertical=max(vertical,eps); end
        draw_reconstruction_curve(horizontal,vertical,method,design);
        if ~isempty(tolerance_field) && options{method}.(tolerance_field)>0
            yline(options{method}.(tolerance_field),':','Color',design.MethodColors(method,:), ...
                'HandleVisibility','off');
        end
    end
    if any(strcmp(quantity,{'residual','residual-time'})) && noise_norm>0
        yline(noise_norm,':k','Actual noise norm','HandleVisibility','off');
    end
    if is_log, set(gca,'YScale','log'); end
    grid on; ylabel(label); title(label); xlabel('Completed step');
    if strcmp(quantity,'residual-time'), xlabel('Solver seconds (including setup)'); end
    legend(names(methods),'Location','best');
end
