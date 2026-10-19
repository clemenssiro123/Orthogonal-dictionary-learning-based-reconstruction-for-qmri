function [pdf_paths,zoom_box] = export_reconstruction_rows(truth,reconstructions,initialization,mask,method_names,design,output_folder)
%EXPORT_RECONSTRUCTION_ROWS Compact, single-page PDFs (MATLAB R2021b).
% [PDF_PATHS,ZOOM_BOX] = EXPORT_RECONSTRUCTION_ROWS(TRUTH,RECONSTRUCTIONS, ...
%     INITIALIZATION,MASK,METHOD_NAMES,DESIGN,OUTPUT_FOLDER)
% TRUTH and each reconstruction cell are N-by-M-by-3 physical maps
% [T1(ms),T2(ms),Density], or N-by-M single images. INITIALIZATION has the
% same size or is empty. MASK is N-by-M; names follow reconstruction order.
% DESIGN is a complete display structure. The destination is created when
% missing. PDF_PATHS lists written files; ZOOM_BOX is the inclusive crop
% [first_row last_row first_column last_column] on the reconstructed grid.
%
% Six files compare T1, T2 and Density maps and absolute errors. Each contains
% one row: ground truth, then RECONSTRUCTIONS in METHOD_NAMES order. The demo
% supplies BLIP, LM, one-step DL and nested DL. The truth error is zero.
% Two more files show initialization maps/errors, ordered T1, T2, Density.
% A single-image input instead writes image_map_comparison.pdf and
% image_error_comparison.pdf. INITIALIZATION=[] omits initialization PDFs.
% All panels share one zoom box. Comparison rows share one colorbar;
% initialization panels have separate colorbars for their different units.
%
% Uses the same DESIGN fields as export_reconstruction_comparison. Optional:
%   RowPanelSize    Longest map side in pixels (default 300).
%   RowPanelGap     Horizontal gap in pixels (default 14).
%   RowShowTitles   Show method/parameter headings (default false).
%   ShowRowFigures  Open compact figure windows (default false); also requires
%                   ShowFigures=true. PDFs are always exported either way.
% Output files replace previous exports with the same names in OUTPUT_FOLDER.
% Returns the exported paths and the actual inclusive ZoomBox coordinates.
%
% Replot an existing experiment without rerunning the reconstruction:
%   s=load(fullfile(run_folder,'results.mat'));
%   s.design.ShowFigures=false;
%   export_reconstruction_rows(s.truth,s.reconstructions,s.u0,s.mask, ...
%       s.method_names,s.design,run_folder);

    %% Validate before creating any figures or files
    [n,m,channels]=size(truth);
    assert(any(channels==[1 3]) && n>=2 && m>=2,'Truth must contain one image or three parameter maps.');
    assert(isequal(size(mask),[n m]),'Mask and parameter maps must have the same size.');
    assert(isempty(initialization) || isequal(size(initialization),size(truth)), ...
        'Initialization must be empty or have the same size as truth.');
    assert(numel(reconstructions)==numel(method_names) && ~isempty(reconstructions), ...
        'Provide one name per reconstruction.');
    for k=1:numel(reconstructions)
        assert(isequal(size(reconstructions{k}),size(truth)), ...
            'Each reconstruction must have the same size as truth.');
    end
    if ~isfield(design,'RowPanelSize'), design.RowPanelSize=300; end
    if ~isfield(design,'RowPanelGap'), design.RowPanelGap=14; end
    if ~isfield(design,'RowShowTitles'), design.RowShowTitles=false; end
    if ~isfield(design,'ShowRowFigures'), design.ShowRowFigures=false; end
    validateattributes(design.RowPanelSize,{'numeric'}, ...
        {'real','finite','scalar','>=',100});
    validateattributes(design.RowPanelGap,{'numeric'}, ...
        {'real','finite','scalar','nonnegative'});
    validateattributes(design.RowShowTitles,{'logical','numeric'}, ...
        {'scalar','binary'});
    validateattributes(design.ShowRowFigures,{'logical','numeric'}, ...
        {'scalar','binary'});
    [zoom_box,design]=comparison_zoom_box(mask,design);
    rows=zoom_box(1):zoom_box(2); cols=zoom_box(3):zoom_box(4);
    if ~isfolder(output_folder), mkdir(output_folder); end

    %% One shared-scale row per parameter and display type
    shown=[{truth} reshape(reconstructions,1,[])];
    headings=[{'Ground truth'} reshape(method_names,1,[])];
    headings(strcmp(headings,'LM'))={'Levenberg-Marquardt'};
    headings(strcmp(headings,'One-step DL'))={'DL (one step)'};
    headings(strcmp(headings,'Nested DL'))={'DL (nested)'};
    stems={'T1','T2','Density'};
    if channels==1, stems={'image'}; end
    % Compact reports currently use unlabeled colorbars and initialization panels.
    labels=repmat({''},1,channels);
    pdf_paths=cell(1,2*channels+2*~isempty(initialization));
    for channel=1:channels
        for is_error=0:1
            values=cell(size(shown));
            for column=1:numel(shown)
                values{column}=display_values(shown{column},truth,channel,is_error,mask,design);
            end
            if is_error
                suffix='error'; limits=design.ErrorLimits(channel,:);
                colors=design.ErrorColors;
            else
                suffix='map'; limits=design.MapLimits(channel,:);
                colors=design.MapColors{channel};
            end
            index=2*channel-1+is_error;
            pdf_paths{index}=fullfile(output_folder,sprintf('%s_%s_comparison.pdf',stems{channel},suffix));
            compact_row(values,headings,repmat(limits,numel(shown),1), ...
                repmat({colors},1,numel(shown)),labels(channel),true,rows,cols,design,pdf_paths{index});
        end
    end

    %% Initialization: three parameters in one row, each with its own scale
    if isempty(initialization), return; end
    for is_error=0:1
        values=cell(1,channels);
        for channel=1:channels
            values{channel}=display_values(initialization,truth,channel,is_error,mask,design);
        end
        if is_error
            filename='initialization_errors.pdf';
            limits=design.ErrorLimits; colors=repmat({design.ErrorColors},1,channels);
        else
            filename='initialization_maps.pdf';
            limits=design.MapLimits; colors=design.MapColors;
        end
        pdf_paths{2*channels+1+is_error}=fullfile(output_folder,filename);
        compact_row(values,labels,limits,colors,labels,false,rows,cols,design,pdf_paths{2*channels+1+is_error});
    end
end

function values=display_values(u,truth,channel,is_error,mask,design)
    % Display copies preserve the physical maps and the saved experiment.
    values=u(:,:,channel);
    if is_error, values=abs(values-truth(:,:,channel)); end
    if design.ZeroBackground, values(~mask)=0; end
end

function compact_row(values,headings,limits,colors,labels,shared_bar,rows,cols,design,path)
    %% Size the page around the images, preserving rectangular map aspect ratios
    [n,m]=size(values{1}); columns=numel(values);
    % The longest image side is RowPanelSize; rectangular maps keep their aspect ratio.
    panel_width=design.RowPanelSize*m/max(n,m);
    panel_height=design.RowPanelSize*n/max(n,m);
    margin=max(8,ceil(design.FontSize)); gap=design.RowPanelGap;
    bar_space=max(85,9*design.FontSize); bar_gap=8; bar_width=12;
    heading_space=0;
    if design.RowShowTitles, heading_space=ceil(2.5*design.FontSize); end
    bars=1; if ~shared_bar, bars=columns; end
    width=2*margin+columns*panel_width+(columns-1)*gap+bars*bar_space;
    height=2*margin+panel_height+heading_space;
    show_figures=design.ShowFigures && design.ShowRowFigures;
    visible='off'; if show_figures, visible='on'; end
    [~,name]=fileparts(path);
    fig=figure('Name',name,'Color','w','Visible',visible,'Units','pixels', ...
        'Position',[40 40 width height],'DefaultAxesFontName',design.FontName, ...
        'DefaultAxesFontSize',design.FontSize,'DefaultTextFontName',design.FontName, ...
        'DefaultTextFontSize',design.FontSize);
    % A white frame fixes a small page margin when exportgraphics crops content.
    annotation(fig,'rectangle',[.001 .001 .998 .998],'Color','w','LineWidth',.1);
    for column=1:columns
        left=margin+(column-1)*(panel_width+gap);
        if ~shared_bar, left=left+(column-1)*bar_space; end
        pos=[left margin panel_width panel_height]./[width height width height];
        ax=axes('Parent',fig,'Units','normalized','Position',pos);
        draw_comparison_panel(ax,values{column},rows,cols,limits(column,:),colors{column},design);
        if design.RowShowTitles
            title(ax,headings{column},'FontWeight','normal','Interpreter','tex');
        end
        if ~shared_bar || column==columns
            cb=colorbar(ax); cb.Units='normalized';
            cb.Position=[left+panel_width+bar_gap margin bar_width panel_height]./[width height width height];
            cb.FontSize=design.FontSize;
            % Individual colorbar labels are omitted when the panel has a heading.
            if shared_bar
                cb.Label.String=labels{1}; cb.Label.Interpreter='tex';
            elseif ~design.RowShowTitles
                cb.Label.String=labels{column}; cb.Label.Interpreter='tex';
            end
            ax.Position=pos; % Undo the axes shrink caused by colorbar creation.
        end
    end
    % Respect the row-figure switch when the shared exporter closes figures.
    row_design=design; row_design.ShowFigures=show_figures;
    export_reconstruction_pdf(fig,path,row_design);
end
