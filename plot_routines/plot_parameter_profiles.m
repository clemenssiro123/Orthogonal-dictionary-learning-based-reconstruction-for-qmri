function fig = plot_parameter_profiles(results,design)
%PLOT_PARAMETER_PROFILES Physical profiles through one image row.
% FIG = PLOT_PARAMETER_PROFILES(RESULTS,DESIGN)
% RESULTS is the four-method results structure, MAT filename, or run folder.
% The three axes compare truth and final reconstructions for T1 (ms),
% T2 (ms), and density. DESIGN.ProfileRow is a one-based row index on the
% reconstructed grid, after downsampling; [] selects round(N/2).
% Optional DESIGN overrides saved display fields. ZeroBackground zeroes
% background samples only in the plotted copies, as in the map figures.
% FIG remains editable; no files are written and no solver is run.
% MATLAB R2021b.
%
% Example: fig=plot_parameter_profiles(s,struct('ProfileRow',40));
% See also PLOT_PARAMETER_MAPS, PLOT_PARAMETER_ZOOM.
    if nargin<2, design=[]; end
    [results,design]=comparison_plot_inputs(results,design);
    truth=results.truth; reconstructions=results.reconstructions;
    mask=results.mask; method_names=results.method_names; map_names=design.MapNames;
    [n,m,~]=size(truth);
    row=design.ProfileRow; if isempty(row), row=round(n/2); end
    validateattributes(row,{'numeric'},{'scalar','integer','>=',1,'<=',n});
    fig=comparison_report_figure('Parameter profiles',design); tiledlayout(1,3,'TileSpacing','compact');
    for j=1:3
        nexttile; hold on; target=truth(:,:,j);
        if design.ZeroBackground, target(~mask)=0; end
        plot(1:m,target(row,:),'k:','LineWidth',2);
        for method=1:4
            values=reconstructions{method}(:,:,j);
            if design.ZeroBackground, values(~mask)=0; end
            draw_reconstruction_curve(1:m,values(row,:),method,design);
        end
        grid on; xlabel('Column'); ylabel(map_names{j}); title(sprintf('Row %d',row));
        if j==1, legend([{'Truth'} method_names],'Location','best'); end
    end
    sgtitle('Physical parameter profiles (display background follows design.ZeroBackground)');

end
