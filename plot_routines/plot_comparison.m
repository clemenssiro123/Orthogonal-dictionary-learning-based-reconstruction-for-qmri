function [fig,zoom_box] = plot_comparison(results,design)
%PLOT_COMPARISON Map/error comparison with linked zoom insets.
% [FIG,ZOOM_BOX] = PLOT_COMPARISON(RESULTS,DESIGN)
% RESULTS is a loaded results structure, MAT filename, or saved-run folder.
% It supplies truth, reconstructions, mask, method_names and saved design.
% Maps have size N-by-M-by-3 in physical [T1(ms),T2(ms),Density] units;
% N-by-M single-image maps are also supported by this comparison grid.
% Optional DESIGN overrides selected saved display fields. ShowFigures
% controls visibility; ZeroBackground masks display copies only.
% FIG is an editable figure. ZOOM_BOX is the actual inclusive crop
% [first_row last_row first_column last_column] on the reconstruction grid.
% No file is written and no solver is run. MATLAB R2021b.
%
% Example: fig=plot_comparison(s,struct('ShowFigures',true));
% See also PLOT_PARAMETER_MAPS, EXPORT_RECONSTRUCTION_COMPARISON.
    if nargin<2, design=[]; end
    [results,design]=comparison_plot_inputs(results,design);
    truth=results.truth; reconstructions=results.reconstructions; mask=results.mask;
    method_names=results.method_names;
    %% 1. Fill optional display settings and check compatible input sizes
    % Only fill missing options; explicitly supplied settings take precedence.
    if ~isfield(design,'ComparisonFigureSize'), design.ComparisonFigureSize=[1100 1320]; end % Use the default portrait figure size when an older design lacks this field.
    [n,m,channels]=size(truth); % Read the reference image height, width and number of parameter channels.
    % assert stops with an explanatory message if a required condition fails.
    % These checks prevent misleading comparisons of differently sized maps.
    assert(any(channels==[1 3]) && n>=2 && m>=2,'Truth must contain one image or three parameter maps.');
    assert(isequal(size(mask),[n m]),'Mask and parameter maps must have the same size.'); % Require a mask with the same two-dimensional image size.
    assert(numel(reconstructions)==numel(method_names) && ~isempty(reconstructions), ... % Require nonempty reconstruction data and one method name per result.
        'Provide one name per reconstruction.'); % Explain a mismatch between results and their names.
    for k=1:numel(reconstructions) % Check the dimensions of every supplied reconstruction.
        assert(isequal(size(reconstructions{k}),size(truth)), ... % Require the same array size as the reference maps.
            'Each reconstruction must have the same size as truth.'); % Explain why differently sized results cannot be plotted together.
    end
    % Require positive figure sizes; shared helpers validate the zoom settings.
    validateattributes(design.ComparisonFigureSize,{'numeric'}, ... % Check the numeric portrait figure-size setting.
        {'real','finite','vector','numel',2,'positive'}); % Require two finite positive values for width and height.
    %% 2. Select one common zoom region for every map and error panel
    [zoom_box,design]=comparison_zoom_box(mask,design); % Reuse the same zoom selection and validation as the compact PDF exporter.
    % Turn the four bounding indices into row/column vectors for direct indexing.
    rows=zoom_box(1):zoom_box(2); cols=zoom_box(3):zoom_box(4); % Turn the inclusive crop boundaries into index vectors for the shared inset.

    %% 3. Arrange the column headings and create the comparison figure
    % Prepend the reference to the methods; its own error panel is therefore zero.
    shown=[{truth} reshape(reconstructions,1,[])]; % Place reference maps first, followed by results in the caller's method order.
    headings=[{'Ground truth'} reshape(method_names,1,[])]; % Construct headings in the same order as the displayed maps.
    % Short publication headings for the demo's two dictionary methods.
    headings(strcmp(headings,'One-step DL'))={'DL (one step)'}; % Shorten the one-step dictionary-learning heading for the portrait grid.
    headings(strcmp(headings,'Nested DL'))={'DL (nested)'}; % Shorten the nested dictionary-learning heading for the portrait grid.
    % Odd rows show maps; the following even rows show absolute errors.
    labels={'','', ... % Define the T1 map/error row labels in milliseconds.
        '','', ... % Define the T2 map/error row labels in milliseconds.
        '',''}; % Define density map/error row labels using the TeX rho symbol.
    % Hidden figures can still be exported for unattended demo runs.
    visible='off'; if design.ShowFigures, visible='on'; end % Choose whether to display the figure window during export.
    fig=figure('Name','Map and error comparison with zoom insets','Color','w', ... % Create a white portrait comparison figure with a descriptive window title.
        'Visible',visible,'Units','pixels','Position',[40 40 design.ComparisonFigureSize], ... % Apply visibility and the configured pixel position/size.
        'DefaultAxesFontName',design.FontName,'DefaultAxesFontSize',design.FontSize, ... % Set the default axes font family and size.
        'DefaultTextFontName',design.FontName,'DefaultTextFontSize',design.FontSize); % Set the default text font family and size.
    % exportgraphics crops to visible content; this white frame preserves
    % page margins around titles, rotated labels and the footer.
    annotation(fig,'rectangle',[.001 .001 .998 .998],'Color','w','LineWidth',.1); % Preserve page margins when exportgraphics crops to the figure's contents.
    % Explicit positions keep ALL panels equally sized, including the ones
    % beside the shared colorbars. Compute in pixels to preserve the layout.
    %% 4. Compute a uniform six-row layout with space for labels and colorbars
    position=get(fig,'Position'); width=position(3); height=position(4); % Read the actual figure dimensions used for layout calculations.
    columns=numel(shown); row_count=2*channels; gap_x=18; gap_y=16; group_gap=12;
    left=65; right=72; top=49; bottom=53; % Reserve outer space for row labels, headings, colorbars and the footer.
    % Fit square panel positions into both available width and available height.
    % The extra group gaps separate the T1, T2 and density row pairs.
    side=min((width-left-right-(columns-1)*gap_x)/columns, ... % Calculate the largest square panel side allowed by the available page width.
        (height-top-bottom-(row_count-1)*gap_y-(channels-1)*group_gap)/row_count);
    assert(side>0,'ComparisonFigureSize is too small for the comparison grid.'); % Reject a figure size that leaves no room for the map grid.
    % Center the complete grid horizontally within the reserved margins.
    left=left+(width-left-right-columns*side-(columns-1)*gap_x)/2; % Center the complete grid horizontally inside the reserved margins.

    %% 5. Draw maps, absolute errors and identical zoom windows
    for row=1:row_count % Alternate map and absolute-error rows for each parameter.
        % Rows 1/2 use channel 1, rows 3/4 channel 2, and rows 5/6 channel 3.
        channel=ceil(row/2); is_error=mod(row,2)==0; % Map each pair of rows to one parameter and mark every even row as an error row.
        % Select these once per row so all methods have identical color scaling.
        if is_error % Select shared display settings for an absolute-error row.
            limits=design.ErrorLimits(channel,:); colors=design.ErrorColors; % Use this parameter's error limits and the error colormap.
        else % Otherwise select settings for a parameter-map row.
            limits=design.MapLimits(channel,:); colors=design.MapColors{channel}; % Use this parameter's map limits and colormap.
        end
        % Figure coordinates start at the bottom; rows are filled from the top.
        bottom_y=height-top-row*side-(row-1)*gap_y-floor((row-1)/2)*group_gap; % Place this row from the top down, including the extra gaps between parameter pairs.
        for column=1:columns % Draw ground truth and all reconstruction columns in order.
            % Convert the calculated pixel position into normalized figure units.
            pos=[left+(column-1)*(side+gap_x) bottom_y side side]./[width height width height]; % Convert this panel's pixel rectangle to normalized figure coordinates.
            ax=axes('Parent',fig,'Units','normalized','Position',pos); % Create an axes with the same panel size as all other columns.
            % Work on a local display array; the stored reconstruction is unchanged.
            values=shown{column}(:,:,channel); % Extract the current parameter from this column's maps.
            % Error is in the same physical units as the corresponding parameter.
            if is_error, values=abs(values-truth(:,:,channel)); end % For error rows, compute the absolute difference from the reference parameter.
            if design.ZeroBackground, values(~mask)=0; end % Optionally mask the background in the display copy only.
            % Draw the full image, crop box, connector and magnified crop together.
            draw_comparison_panel(ax,values,rows,cols,limits,colors,design); % Reuse the common full-map/inset drawing helper.
            % A single title at the top identifies each complete method column.
            if row==1 % Put a column heading above the first row only.
                title(ax,headings{column},'FontWeight','normal', ... % Draw the current method name without bold text.
                    'FontSize',design.FontSize+3,'Interpreter','none'); % Use a slightly larger heading font and keep the name as literal text.
            end
            % A single rotated label at the left identifies each complete row.
            if column==1 % Put the parameter/error label beside the first column only.
                text(ax,-.10,.5,labels{row},'Units','normalized','Rotation',90, ... % Place the label left of the image and rotate it vertically.
                    'HorizontalAlignment','center','VerticalAlignment','bottom', ... % Center the label along the row while anchoring its edge beside the panel.
                    'FontSize',design.FontSize+1,'Interpreter','tex','Clipping','off'); % Allow TeX subscripts and let the label extend outside the image axes.
            end
            % Share the last panel's colorbar across its entire row.
            if column==columns % Place one shared colorbar beside the last panel in each row.
                cb=colorbar(ax); % Create a bar linked to this row's color limits and colormap.
                cb.Units='normalized'; % Use normalized figure units for manual colorbar placement.
                cb.Position=[pos(1)+pos(3)+9/width pos(2) 10/width pos(4)]; % Place the bar beside the last panel and match the panel height.
                cb.FontSize=design.FontSize-1; % Use slightly smaller tick text to keep the portrait grid compact.
                % Creating a colorbar can shrink axes; restore the common panel size.
                ax.Position=pos; % Restore panel geometry after MATLAB makes room for the colorbar.
            end
        end
    end
    %% 6. Record the crop coordinates
    % The footer records both the crop and whether the body mask was displayed.
    background='Background shown'; % Prepare the footer statement for an unmasked display.
    if design.ZeroBackground, background='Background masked for display'; end % Replace that statement when the displayed background has been zeroed.
    caption=sprintf(['Zoom: rows %d-%d, columns %d-%d (reconstruction grid).  %s.\n' ... % Build the footer with crop indices and the background-display description.
        'Absolute errors in physical units; shared color limits within each row; values above the limits saturate.'], ... % Explain error units, common scales and saturation above the displayed limits.
        zoom_box,background); % Insert the selected crop coordinates and background description into the caption.
    annotation(fig,'textbox',[.045 .006 .94 .03],'String',caption, ... % Place the footer in a textbox along the bottom of the page.
        'EdgeColor','none','FontName',design.FontName,'FontSize',design.FontSize-1, ... % Remove its border and use the configured font with a slightly smaller size.
        'Interpreter','none','VerticalAlignment','middle','HorizontalAlignment','center'); % Render literal text and center the footer within its textbox.
end
