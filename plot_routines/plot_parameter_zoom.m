function fig = plot_parameter_zoom(results,design)
%PLOT_PARAMETER_ZOOM Enlarged maps for an explicit or automatic crop.
% FIG = PLOT_PARAMETER_ZOOM(RESULTS,DESIGN)
% RESULTS is the four-method results structure, MAT filename, or run folder.
% Optional DESIGN overrides saved fields. ZoomBox is the inclusive vector
% [first_row last_row first_column last_column] on the reconstructed grid.
% The crop must span at least two rows and two columns. [] uses the same
% automatic body-region selection as PLOT_COMPARISON.
% Shows truth and all four reconstructions in physical units. Color limits
% and ZeroBackground affect display copies only.
% FIG remains editable; no files are written and no solver is run.
% MATLAB R2021b.
% See also PLOT_COMPARISON, PLOT_PARAMETER_PROFILES.
    if nargin<2, design=[]; end
    [results,design]=comparison_plot_inputs(results,design);
    truth=results.truth; reconstructions=results.reconstructions;
    mask=results.mask; method_names=results.method_names; map_names=design.MapNames;
    [design.ZoomBox,design]=comparison_zoom_box(mask,design);
    box=design.ZoomBox; rows=box(1):box(2); columns=box(3):box(4);
    fig=comparison_report_figure('Enlarged parameter maps',design); tiledlayout(3,5,'TileSpacing','compact');
    shown=[{truth} reconstructions]; headings=[{'Truth'} method_names];
    for j=1:3
        for column=1:5
            ax=nexttile;
            show_reconstruction_map(ax,shown{column}(rows,columns,j),mask(rows,columns), ...
                design.MapLimits(j,:),design.MapColors{j},design);
            title([headings{column} ' | ' map_names{j}]); if column==5, drawnow; colorbar(ax); end
        end
    end
    sgtitle(sprintf('Zoom: rows %d-%d, columns %d-%d',box));
end
