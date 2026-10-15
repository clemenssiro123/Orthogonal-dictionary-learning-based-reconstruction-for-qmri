function zoom_box = export_reconstruction_comparison(truth,reconstructions,mask,method_names,design,path,append) % Write or append the portrait comparison page to path and return the chosen crop.
%EXPORT_RECONSTRUCTION_COMPARISON Map/error grid with linked zoom insets.
% MATLAB R2021b. Display only: no interpolation toolbox or solver rerun needed.
% Columns are ground truth followed by the reconstruction methods. Each
% parameter has a map row and an absolute-error row, with one colorbar per row.
%
% Uses the demo's MapColors, MapLimits, ErrorColors, ErrorLimits, FontName,
% FontSize, ShowFigures, ZeroBackground and Resolution fields. Optional:
%   ZoomBox               [first_row last_row first_col last_col], inclusive,
%                         on the reduced reconstruction grid. [] selects a
%                         lower-left region within the body's bounding box.
%   ZoomInsetFraction     Maximum inset side / map side (default .43).
%   ZoomColor             ROI, connector and inset border RGB (default red).
%   ComparisonFigureSize  Portrait figure size in pixels (default [1100 1320]).
% APPEND defaults to false; true adds the comparison to an existing PDF.
% Returns the actual ZoomBox used so automatic selections can be reproduced.
% TRUTH and each cell of RECONSTRUCTIONS have size N-by-M-by-3, with
% channels [T1,T2,rho] in physical units. MASK is an N-by-M logical body mask.
% DESIGN controls presentation; PATH is the destination PDF filename.
% All columns use the same crop and each row uses the same color scale,
% so differences between methods remain directly comparable.
%
% Replot a saved run without recomputing any reconstruction:
%   s=load(fullfile(run_folder,'results.mat'));
%   s.design.ShowFigures=false;
%   export_reconstruction_comparison(s.truth,s.reconstructions,s.mask, ...
%       s.method_names,s.design,fullfile(run_folder,'comparison_with_zooms.pdf'));

    %% 1. Fill optional display settings and check compatible input sizes
    % nargin counts supplied inputs; omitting APPEND starts a new PDF.
    if nargin<7, append=false; end % Start a new PDF by default when the optional append argument is omitted.
    % Only fill missing options; explicitly supplied settings take precedence.
    if ~isfield(design,'ComparisonFigureSize'), design.ComparisonFigureSize=[1100 1320]; end % Use the default portrait figure size when an older design lacks this field.
    [n,m,channels]=size(truth); % Read the reference image height, width and number of parameter channels.
    % assert stops with an explanatory message if a required condition fails.
    % These checks prevent misleading comparisons of differently sized maps.
    assert(channels==3 && n>=2 && m>=2,'Truth must contain three 2-D parameter maps.'); % Require three parameter maps with at least two rows and columns.
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
    labels={'T_1 maps (ms)','T_1 absolute error (ms)', ... % Define the T1 map/error row labels in milliseconds.
        'T_2 maps (ms)','T_2 absolute error (ms)', ... % Define the T2 map/error row labels in milliseconds.
        '\rho maps','\rho absolute error'}; % Define density map/error row labels using the TeX rho symbol.
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
    columns=numel(shown); gap_x=18; gap_y=16; group_gap=12; % Count columns and set pixel gaps between columns, rows and parameter groups.
    left=65; right=72; top=49; bottom=53; % Reserve outer space for row labels, headings, colorbars and the footer.
    % Fit square panel positions into both available width and available height.
    % The extra group gaps separate the T1, T2 and density row pairs.
    side=min((width-left-right-(columns-1)*gap_x)/columns, ... % Calculate the largest square panel side allowed by the available page width.
        (height-top-bottom-5*gap_y-2*group_gap)/6); % Also constrain the side to fit all six rows and their vertical gaps.
    assert(side>0,'ComparisonFigureSize is too small for the comparison grid.'); % Reject a figure size that leaves no room for the map grid.
    % Center the complete grid horizontally within the reserved margins.
    left=left+(width-left-right-columns*side-(columns-1)*gap_x)/2; % Center the complete grid horizontally inside the reserved margins.

    %% 5. Draw maps, absolute errors and identical zoom windows
    for row=1:6 % Draw the six rows: T1 map/error, T2 map/error, and density map/error.
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
    %% 6. Record the crop coordinates and write or append the PDF page
    % The footer records both the crop and whether the body mask was displayed.
    background='Background shown'; % Prepare the footer statement for an unmasked display.
    if design.ZeroBackground, background='Background masked for display'; end % Replace that statement when the displayed background has been zeroed.
    caption=sprintf(['Zoom: rows %d-%d, columns %d-%d (reconstruction grid).  %s.\n' ... % Build the footer with crop indices and the background-display description.
        'Absolute errors in physical units; shared color limits within each row; values above the limits saturate.'], ... % Explain error units, common scales and saturation above the displayed limits.
        zoom_box,background); % Insert the selected crop coordinates and background description into the caption.
    annotation(fig,'textbox',[.045 .006 .94 .03],'String',caption, ... % Place the footer in a textbox along the bottom of the page.
        'EdgeColor','none','FontName',design.FontName,'FontSize',design.FontSize-1, ... % Remove its border and use the configured font with a slightly smaller size.
        'Interpreter','none','VerticalAlignment','middle','HorizontalAlignment','center'); % Render literal text and center the footer within its textbox.
    % Finish graphics layout before exporting; APPEND preserves earlier PDF pages.
    drawnow; % Wait for the graphics layout to finish before exporting.
    exportgraphics(fig,path,'ContentType','vector','BackgroundColor','white', ... % WRITE THE PDF HERE: export the portrait figure to the caller's path.
        'Resolution',design.Resolution,'Append',append); % Use the configured resolution and either append a page or replace the PDF.
    % Hidden figures have served their purpose once the PDF has been written.
    if ~design.ShowFigures, close(fig); end % Close hidden figures after export; retain visible figures for inspection.
end
