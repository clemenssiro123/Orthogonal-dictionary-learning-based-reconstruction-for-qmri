function fig = plot_parameter_errors(results,design)
%PLOT_PARAMETER_ERRORS Absolute errors with shared color limits.
% FIG = PLOT_PARAMETER_ERRORS(RESULTS,DESIGN) recreates errors.pdf.
% RESULTS is the four-method results structure, MAT filename, or run folder.
% Images show abs(reconstruction-truth), in ms for T1/T2 and original density
% units. Titles use the saved relative L2 body errors in RESULTS.comparison;
% these percentages differ from the absolute errors shown in the images.
% Optional DESIGN overrides saved fields such as ErrorLimits and
% ZeroBackground. Saturated colors do not clip the underlying error values.
% FIG remains editable; no files are written and no solver is run.
% MATLAB R2021b.
% See also PLOT_COMPARISON, PLOT_CONVERGENCE.
    if nargin<2, design=[]; end
    [results,design]=comparison_plot_inputs(results,design);
    truth=results.truth; reconstructions=results.reconstructions;
    mask=results.mask; method_names=results.method_names; map_names=design.MapNames;
    relative_error=results.comparison{:,{'T1Error_percent','T2Error_percent','DensityError_percent'}};
    fig=comparison_report_figure('Absolute reconstruction errors',design);
    % Absolute errors are in ms for T1/T2 and the density map's original units.
    tiledlayout(3,4,'TileSpacing','compact','Padding','compact');
    for j=1:3
        for method=1:4
            ax=nexttile;
            show_reconstruction_map(ax,abs(reconstructions{method}(:,:,j)-truth(:,:,j)),mask, ...
                design.ErrorLimits(j,:),design.ErrorColors,design);
            title(sprintf('%s | %s\nBody error %.2f%%',method_names{method},map_names{j},relative_error(method,j)));
            if method==4, drawnow; colorbar(ax); end
        end
    end
    sgtitle('Absolute errors in physical units; common limits per row (larger errors saturate)');

end
