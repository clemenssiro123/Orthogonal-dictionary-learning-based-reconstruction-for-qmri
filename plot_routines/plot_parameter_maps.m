function fig = plot_parameter_maps(results,design)
%PLOT_PARAMETER_MAPS Truth, common start, and reconstructed maps.
% FIG = PLOT_PARAMETER_MAPS(RESULTS,DESIGN) recreates the first maps.pdf page.
% RESULTS is the four-method results structure, MAT filename, or run folder.
% Columns contain truth, coarse-BLIP start, BLIP, LM, one-step DL and nested
% DL. Rows contain T1 (ms), T2 (ms), and density in its original units.
% The acquisition caption uses RESULTS.experiment, not current presets.
% Optional DESIGN overrides saved fields, for example ShowFigures or
% MapLimits. Color limits and background masking affect display copies only.
% FIG remains editable; no files are written and no solver is run.
% MATLAB R2021b.
% See also PLOT_COMPARISON, EXPORT_COMPARISON_REPORTS.
    if nargin<2, design=[]; end
    [results,design]=comparison_plot_inputs(results,design);
    truth=results.truth; u0=results.u0; reconstructions=results.reconstructions;
    mask=results.mask; method_names=results.method_names; map_names=design.MapNames;
    [n,m,~]=size(truth); L=results.experiment.Frames;
    r=results.experiment.Undersampling; noise_std=results.experiment.NoiseStd;
    fig=comparison_report_figure('Parameter maps',design); tiledlayout(3,6,'TileSpacing','compact','Padding','compact');
    shown=[{truth,u0} reconstructions]; headings=[{'Truth','Coarse BLIP start'} method_names];
    for j=1:3
        for column=1:6
            ax=nexttile;
            show_reconstruction_map(ax,shown{column}(:,:,j),mask,design.MapLimits(j,:),design.MapColors{j},design);
            title(sprintf('%s\n%s',headings{column},map_names{j}));
            if column==6, drawnow; colorbar(ax); end
        end
    end
    sgtitle(sprintf('Shared start | %dx%d, L=%d, r=%d, noise SD=%.2g per component',n,m,L,r,noise_std));

end
