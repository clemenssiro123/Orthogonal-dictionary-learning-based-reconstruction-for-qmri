function [pdf_paths,zoom_box] = export_reconstruction_rows(truth,reconstructions,initialization,mask,method_names,design,output_folder) % Write the eight compact PDFs to output_folder and return their paths and shared zoom box.
%EXPORT_RECONSTRUCTION_ROWS Eight compact, single-page PDFs (MATLAB R2021b).
% Six files compare T1, T2 and Density maps and absolute errors. Each contains
% one row: ground truth, then RECONSTRUCTIONS in METHOD_NAMES order. The demo
% supplies BLIP, LM, one-step DL and nested DL. The truth error is zero.
% Two more files show initialization maps/errors, ordered T1, T2, Density.
% All panels share one zoom box. Comparison rows share one colorbar;
% initialization panels have separate colorbars for their different units.
%
% Uses the same DESIGN fields as export_reconstruction_comparison. Optional:
%   RowPanelSize    Longest map side in pixels (default 300).
%   RowPanelGap     Horizontal gap in pixels (default 14).
%   RowShowTitles   Show method/parameter headings (default true).
%   ShowRowFigures  Open compact figure windows (default false); also requires
%                   ShowFigures=true. PDFs are always exported either way.
% Output files replace previous exports with the same names in OUTPUT_FOLDER.
% Returns the eight paths and the actual inclusive ZoomBox coordinates.
%
% Replot an existing experiment without rerunning the reconstruction:
%   s=load(fullfile(run_folder,'results.mat'));
%   s.design.ShowFigures=false;
%   export_reconstruction_rows(s.truth,s.reconstructions,s.u0,s.mask, ...
%       s.method_names,s.design,run_folder);

    %% Validate before creating any figures or files
    [n,m,channels]=size(truth); % Read image height, width and parameter-channel count from the reference maps.
    assert(channels==3 && n>=2 && m>=2,'Truth must contain three 2-D parameter maps.'); % Require T1, T2 and density maps with at least two rows and columns.
    assert(isequal(size(mask),[n m]),'Mask and parameter maps must have the same size.'); % Require one mask entry per image pixel.
    assert(isequal(size(initialization),size(truth)), ... % Check that initialization and reference arrays have identical dimensions.
        'Initialization must have the same size as truth.'); % Report incompatible initialization dimensions before writing files.
    assert(numel(reconstructions)==numel(method_names) && ~isempty(reconstructions), ... % Require at least one reconstruction and one method name for each result.
        'Provide one name per reconstruction.'); % Report an empty reconstruction list or a mismatch in method-name count.
    for k=1:numel(reconstructions) % Validate each reconstruction array in turn.
        assert(isequal(size(reconstructions{k}),size(truth)), ... % Compare this result's dimensions with the reference array.
            'Each reconstruction must have the same size as truth.'); % Report any result that cannot be compared pixel by pixel.
    end
    if ~isfield(design,'RowPanelSize'), design.RowPanelSize=300; end % Default the longest displayed map side to 300 pixels for older saved designs.
    if ~isfield(design,'RowPanelGap'), design.RowPanelGap=14; end % Default the horizontal gap between maps to 14 pixels.
    if ~isfield(design,'RowShowTitles'), design.RowShowTitles=true; end % Show method or parameter headings unless explicitly disabled.
    if ~isfield(design,'ShowRowFigures'), design.ShowRowFigures=false; end % Export compact PDFs without opening windows by default, including for older saved runs.
    validateattributes(design.RowPanelSize,{'numeric'}, ... % Check the numeric map-size setting.
        {'real','finite','scalar','>=',100}); % Require a finite scalar size of at least 100 pixels.
    validateattributes(design.RowPanelGap,{'numeric'}, ... % Check the numeric inter-panel spacing setting.
        {'real','finite','scalar','nonnegative'}); % Require a finite nonnegative scalar gap in pixels.
    validateattributes(design.RowShowTitles,{'logical','numeric'}, ... % Check the logical or numeric heading switch.
        {'scalar','binary'}); % Accept only a scalar false/true or zero/one value.
    validateattributes(design.ShowRowFigures,{'logical','numeric'}, ... % Check the compact-figure visibility switch before creating files.
        {'scalar','binary'}); % Accept only a scalar false/true or zero/one value.
    [zoom_box,design]=comparison_zoom_box(mask,design); % Choose and validate one common crop and fill missing inset settings.
    rows=zoom_box(1):zoom_box(2); cols=zoom_box(3):zoom_box(4); % Expand the inclusive crop bounds into row and column index vectors.
    if ~isfolder(output_folder), mkdir(output_folder); end % Create the destination directory if this is its first export.

    %% One shared-scale row per parameter and display type
    shown=[{truth} reshape(reconstructions,1,[])]; % Put ground truth first and retain the supplied reconstruction order.
    headings=[{'Ground truth'} reshape(method_names,1,[])]; % Build the corresponding left-to-right column headings.
    headings(strcmp(headings,'LM'))={'Levenberg-Marquardt'}; % Expand the LM abbreviation for the displayed heading.
    headings(strcmp(headings,'One-step DL'))={'DL (one step)'}; % Use a short, consistent heading for one-step dictionary learning.
    headings(strcmp(headings,'Nested DL'))={'DL (nested)'}; % Use a short, consistent heading for nested dictionary learning.
    stems={'T1','T2','Density'}; % Define the three parameter prefixes used in the PDF filenames.
    map_labels={'T_1 (ms)','T_2 (ms)','Density'}; % Define map colorbar labels; TeX underscores produce subscripts.
    error_labels={'T_1 absolute error (ms)','T_2 absolute error (ms)','Density absolute error'}; % Define error labels in the same physical units as the corresponding maps.
    pdf_paths=cell(1,8); % Reserve one returned filename for each of the eight PDFs.
    for channel=1:3 % Export T1, T2 and density in channel order.
        for is_error=0:1 % For each channel, export the map first and its absolute error second.
            values=cell(size(shown)); % Allocate a display array for every comparison column.
            for column=1:numel(shown) % Prepare ground truth and each reconstruction in the same order as the headings.
                values{column}=display_values(shown{column},truth,channel,is_error,mask,design); % Select the channel, optionally compute absolute error, and apply display masking.
            end
            if is_error % Choose settings for the absolute-error PDF.
                suffix='error'; limits=design.ErrorLimits(channel,:); % Select the error filename suffix and this channel's shared error limits.
                colors=design.ErrorColors; label=error_labels{channel}; % Select the common error colormap and the physical-unit error label.
            else % Otherwise choose settings for the parameter-map PDF.
                suffix='map'; limits=design.MapLimits(channel,:); % Select the map filename suffix and this channel's shared map limits.
                colors=design.MapColors{channel}; label=map_labels{channel}; % Select this parameter's map colormap and physical-unit label.
            end
            index=2*channel-1+is_error; % Assign slots 1/2 to T1, 3/4 to T2, and 5/6 to density.
            pdf_paths{index}=fullfile(output_folder,sprintf('%s_%s_comparison.pdf',stems{channel},suffix)); % Build the full output filename, for example output_folder/T1_map_comparison.pdf.
            compact_row(values,headings,repmat(limits,numel(shown),1), ... % Draw the row with the same numeric color limits repeated for every method.
                repmat({colors},1,numel(shown)),{label},true,rows,cols,design,pdf_paths{index}); % Use one shared colorbar and export this page to the computed PDF path.
        end
    end

    %% Initialization: three parameters in one row, each with its own scale
    for is_error=0:1 % Export initialization maps first, then initialization absolute errors.
        values=cell(1,3); % Allocate one display panel for each of the three parameters.
        for channel=1:3 % Prepare T1, T2 and density in left-to-right order.
            values{channel}=display_values(initialization,truth,channel,is_error,mask,design); % Extract the initial map or its absolute error and apply display-only masking.
        end
        if is_error % Choose settings for the initialization error page.
            filename='initialization_errors.pdf'; labels=error_labels; % Select the initialization error filename and its three parameter headings.
            limits=design.ErrorLimits; colors=repmat({design.ErrorColors},1,3); % Use each parameter's error limits and the common error colormap.
        else % Otherwise choose settings for the initialization map page.
            filename='initialization_maps.pdf'; labels=map_labels; % Select the initialization map filename and its three parameter headings.
            limits=design.MapLimits; colors=design.MapColors; % Use each parameter's map limits and its own colormap.
        end
        pdf_paths{7+is_error}=fullfile(output_folder,filename); % Store initialization map/error paths in the final two return slots.
        compact_row(values,labels,limits,colors,labels,false,rows,cols,design,pdf_paths{7+is_error}); % Export a three-panel page with a separate colorbar for each parameter.
    end
end

function values=display_values(u,truth,channel,is_error,mask,design) % Prepare one display array without modifying the input maps.
    % Display copies preserve the physical maps and the saved experiment.
    values=u(:,:,channel); % Copy the requested T1, T2 or density channel.
    if is_error, values=abs(values-truth(:,:,channel)); end % For error panels, compute the absolute pixelwise difference from ground truth.
    if design.ZeroBackground, values(~mask)=0; end % Optionally zero pixels outside the body mask in this display copy only.
end

function compact_row(values,headings,limits,colors,labels,shared_bar,rows,cols,design,path) % Arrange one narrow figure and write its single-page PDF to path.
    %% Size the page around the images, preserving rectangular map aspect ratios
    [n,m]=size(values{1}); columns=numel(values); % Read the common map shape and the number of columns to lay out.
    panel_width=design.RowPanelSize*m/max(n,m); % Scale map width so the longer image dimension equals RowPanelSize.
    panel_height=design.RowPanelSize*n/max(n,m); % Scale map height by the same factor to preserve anatomy proportions.
    margin=max(8,ceil(design.FontSize)); gap=design.RowPanelGap; % Reserve a small outer margin and read the configured inter-panel gap.
    bar_space=max(85,9*design.FontSize); bar_gap=8; bar_width=12; % Reserve colorbar/label space; set the map-to-bar gap and bar width in pixels.
    heading_space=0; % Start without extra vertical space for headings.
    if design.RowShowTitles, heading_space=ceil(2.5*design.FontSize); end % Reserve heading space proportional to font size only when headings are enabled.
    bars=1; if ~shared_bar, bars=columns; end % Reserve one shared colorbar, or one colorbar per parameter for initialization.
    width=2*margin+columns*panel_width+(columns-1)*gap+bars*bar_space; % Add margins, map widths, gaps and colorbar spaces to obtain page width.
    height=2*margin+panel_height+heading_space; % Add margins, image height and heading space to obtain page height.
    show_figures=design.ShowFigures && design.ShowRowFigures; % Open compact windows only when both the master switch and the compact-figure switch are enabled.
    visible='off'; if show_figures, visible='on'; end % Export using hidden figures by default while retaining the option to inspect them on screen.
    [~,name]=fileparts(path); % Use the PDF filename without its extension as the figure-window name.
    fig=figure('Name',name,'Color','w','Visible',visible,'Units','pixels', ... % Create a white figure with pixel-based dimensions and the chosen visibility.
        'Position',[40 40 width height],'DefaultAxesFontName',design.FontName, ... % Set its screen position, computed page size and axes font family.
        'DefaultAxesFontSize',design.FontSize,'DefaultTextFontName',design.FontName, ... % Apply the configured axes font size and text font family.
        'DefaultTextFontSize',design.FontSize); % Apply the configured text font size to new text objects.
    if ~show_figures, cleanup=onCleanup(@() close(fig)); end % Automatically close a hidden figure on return or error; retain explicitly requested visible figures.
    % A white frame fixes a small page margin when exportgraphics crops content.
    annotation(fig,'rectangle',[.001 .001 .998 .998],'Color','w','LineWidth',.1); % Keep a nearly full-page white frame so PDF cropping preserves the small margin.
    for column=1:columns % Draw every panel from left to right.
        left=margin+(column-1)*(panel_width+gap); % Compute this panel's horizontal starting point from its index and gap.
        if ~shared_bar, left=left+(column-1)*bar_space; end % For initialization pages, also skip space occupied by earlier colorbars.
        pos=[left margin panel_width panel_height]./[width height width height]; % Convert pixel geometry to normalized [left bottom width height] coordinates.
        ax=axes('Parent',fig,'Units','normalized','Position',pos); % Create an axes at the calculated panel position.
        draw_comparison_panel(ax,values{column},rows,cols,limits(column,:),colors{column},design); % Draw the full map, source-region box, connector and magnified inset.
        if design.RowShowTitles % Add a heading only when requested by the design settings.
            title(ax,headings{column},'FontWeight','normal','Interpreter','tex'); % Center the method/parameter heading above the image; allow TeX subscripts.
        end
        if ~shared_bar || column==columns % Add a bar beside every initialization panel or beside the last comparison panel.
            cb=colorbar(ax); cb.Units='normalized'; % Create the axes colorbar and use normalized figure coordinates for placement.
            cb.Position=[left+panel_width+bar_gap margin bar_width panel_height]./[width height width height]; % Place a thin bar immediately to the right, matching the displayed map height.
            cb.FontSize=design.FontSize; % Match colorbar tick text to the configured font size.
            % Comparison pages identify the quantity beside their shared scale.
            % Initialization headings already identify each scale; add a label
            % only when headings are disabled.
            if shared_bar % Label the shared comparison scale with the quantity and its units.
                cb.Label.String=labels{1}; cb.Label.Interpreter='tex'; % Set the common label and enable TeX subscripts such as T_1.
            elseif ~design.RowShowTitles % For title-free initialization pages, put parameter identification beside each bar.
                cb.Label.String=labels{column}; cb.Label.Interpreter='tex'; % Use this parameter's label with TeX formatting for subscripts.
            end
            ax.Position=pos; % Undo any axes shrinking caused by colorbar creation to preserve equal map sizes.
        end
    end
    drawnow; % Complete pending graphics layout before exporting the figure.
    exportgraphics(fig,path,'ContentType','vector','BackgroundColor','white', ... % WRITE THE PDF HERE: export this figure to path with vector graphics and a white background.
        'Resolution',design.Resolution); % Supply the configured raster resolution where applicable; text/lines remain vector graphics.
end
