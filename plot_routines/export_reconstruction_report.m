function export_reconstruction_report(truth,reconstructions,initialization,mask,method_names,design,output_folder)
%EXPORT_RECONSTRUCTION_REPORT Common maps, errors, compact rows and profiles.
% EXPORT_RECONSTRUCTION_REPORT(TRUTH,RECONSTRUCTIONS,INITIALIZATION,MASK, ...
%     METHOD_NAMES,DESIGN,OUTPUT_FOLDER)
% TRUTH and each reconstruction cell contain N-by-M-by-3 physical maps
% [T1(ms),T2(ms),Density], or N-by-M single images. MASK is N-by-M, and
% METHOD_NAMES follows the reconstruction order. INITIALIZATION has the same
% map size, or [] when no initial maps are available (legacy BLIP).
% DESIGN is a complete structure from reconstruction_plot_design or a saved
% run. It controls colors, scales, display masking, zooms, profiles and fonts.
% OUTPUT_FOLDER is created if needed; PDFs with matching names are replaced.
% ShowFigures=false closes hidden figures after export. Stored maps remain
% unchanged, and scalar errors are computed on MASK pixels before display.
% Convergence diagnostics are exported separately by each standalone demo.
% For the full saved four-method report, use EXPORT_COMPARISON_REPORTS.
% MATLAB R2021b.
% See also EXPORT_COMPARISON_REPORTS, EXPORT_RECONSTRUCTION_ROWS.
    [n,m,channels]=size(truth);
    assert(any(channels==[1 3]) && n>=2 && m>=2,'Provide one image or three parameter maps.');
    assert(isequal(size(mask),[n m]),'Mask and maps must have the same size.');
    reconstructions=reshape(reconstructions,1,[]);
    method_names=reshape(method_names,1,[]);
    assert(~isempty(reconstructions) && numel(reconstructions)==numel(method_names), ...
        'Provide one name per reconstruction.');
    for method=1:numel(reconstructions)
        assert(isequal(size(reconstructions{method}),size(truth)), ...
            'Each reconstruction must have the same size as truth.');
    end
    assert(isempty(initialization) || isequal(size(initialization),size(truth)), ...
        'Initialization must be empty or have the same size as truth.');
    [~,design]=comparison_zoom_box(mask,design);
    row=design.ProfileRow; if isempty(row), row=round(n/2); end
    validateattributes(row,{'numeric'},{'scalar','integer','>=',1,'<=',n});
    if ~isfolder(output_folder), mkdir(output_folder); end
    shown=[{truth} reconstructions]; headings=[{'Truth'} method_names];
    if ~isempty(initialization)
        shown=[{truth,initialization} reconstructions];
        headings=[{'Truth','Initial'} method_names];
    end

    fig=reconstruction_report_figure('Parameter maps',design);
    tiledlayout(channels,numel(shown),'TileSpacing','compact','Padding','compact');
    for j=1:channels
        for column=1:numel(shown)
            ax=nexttile;
            show_map(ax,shown{column}(:,:,j),mask,design.MapLimits(j,:),design.MapColors{j},design);
            title(sprintf('%s\n%s',headings{column},design.MapNames{j}),'Interpreter','none');
            if column==numel(shown), colorbar(ax); end
        end
    end
    sgtitle('Physical maps; common color limits per row');
    export_reconstruction_pdf(fig,fullfile(output_folder,'maps.pdf'),design);
    export_reconstruction_comparison(truth,reconstructions,mask,method_names, ...
        design,fullfile(output_folder,'maps.pdf'),true);
    export_reconstruction_rows(truth,reconstructions,initialization,mask, ...
        method_names,design,output_folder);

    fig=reconstruction_report_figure('Absolute reconstruction errors',design);
    tiledlayout(channels,numel(reconstructions),'TileSpacing','compact','Padding','compact');
    relative_error=zeros(numel(reconstructions),channels);
    for j=1:channels
        target=truth(:,:,j);
        for method=1:numel(reconstructions)
            difference=reconstructions{method}(:,:,j)-target;
            relative_error(method,j)=100*norm(difference(mask))/max(norm(target(mask)),eps);
            ax=nexttile;
            show_map(ax,abs(difference),mask,design.ErrorLimits(j,:),design.ErrorColors,design);
            title(sprintf('%s | %s\nMask error %.2f%%',method_names{method}, ...
                design.MapNames{j},relative_error(method,j)),'Interpreter','none');
            if method==numel(reconstructions), colorbar(ax); end
        end
    end
    sgtitle('Absolute errors in physical units; larger errors saturate');
    export_reconstruction_pdf(fig,fullfile(output_folder,'errors.pdf'),design);

    profiles_pdf=fullfile(output_folder,'profiles_and_summary.pdf');
    fig=reconstruction_report_figure('Parameter profiles',design);
    tiledlayout(1,channels,'TileSpacing','compact','Padding','compact');
    for j=1:channels
        nexttile; hold on; target=truth(:,:,j);
        if design.ZeroBackground, target(~mask)=0; end
        plot(1:m,target(row,:),'k:','LineWidth',2);
        for method=1:numel(reconstructions)
            values=reconstructions{method}(:,:,j);
            if design.ZeroBackground, values(~mask)=0; end
            plot(1:m,values(row,:),'Color',design.MethodColors(method,:), ...
                'LineStyle',design.LineStyles{method},'LineWidth',design.LineWidth);
        end
        grid on; xlabel('Column'); ylabel(design.MapNames{j},'Interpreter','none');
        title(sprintf('Row %d',row));
        if j==1, legend([{'Truth'} method_names],'Location','best','Interpreter','none'); end
    end
    sgtitle('Physical profiles; background follows design.ZeroBackground');
    export_reconstruction_pdf(fig,profiles_pdf,design);
    if ~isempty(design.ZoomBox)
        box=design.ZoomBox; rows=box(1):box(2); cols=box(3):box(4);
        zoomed=[{truth} reconstructions]; zoom_names=[{'Truth'} method_names];
        fig=reconstruction_report_figure('Enlarged parameter maps',design);
        tiledlayout(channels,numel(zoomed),'TileSpacing','compact','Padding','compact');
        for j=1:channels
            for column=1:numel(zoomed)
                ax=nexttile;
                show_map(ax,zoomed{column}(rows,cols,j),mask(rows,cols), ...
                    design.MapLimits(j,:),design.MapColors{j},design);
                title([zoom_names{column} ' | ' design.MapNames{j}],'Interpreter','none');
                if column==numel(zoomed), colorbar(ax); end
            end
        end
        sgtitle(sprintf('Zoom: rows %d-%d, columns %d-%d',box));
        export_reconstruction_pdf(fig,profiles_pdf,design,true);
    end
    fig=reconstruction_report_figure('Reconstruction summary',design);
    ax=axes('Parent',fig,'Position',[.05 .05 .9 .9]); axis(ax,'off');
    lines_out={sprintf('RECONSTRUCTION ERRORS | %dx%d pixels',n,m), ...
        'Relative L2 errors (%) on mask pixels, using the unmodified maps.', ' '};
    for method=1:numel(reconstructions)
        lines_out{end+1}=method_names{method}; %#ok<AGROW>
        for j=1:channels
            difference=reconstructions{method}(:,:,j)-truth(:,:,j);
            lines_out{end+1}=sprintf('  %-12s  error %8.3f%%   RMSE %10.4g', ...
                design.MapNames{j},relative_error(method,j),sqrt(mean(difference(mask).^2))); %#ok<AGROW>
        end
    end
    lines_out=[lines_out {' ','Color limits and background masking affect display copies only.'}];
    text(ax,0,1,lines_out,'VerticalAlignment','top','Interpreter','none', ...
        'FontName',design.ReportFontName,'FontSize',design.ReportFontSize);
    export_reconstruction_pdf(fig,profiles_pdf,design,true);
    fprintf('\nPDF reports saved in:\n%s\n',output_folder);
end

function show_map(ax,values,mask,limits,colors,design)
    if design.ZeroBackground, values(~mask)=0; end
    imagesc(ax,values,limits); axis(ax,'image'); axis(ax,'off'); colormap(ax,colors);
end
