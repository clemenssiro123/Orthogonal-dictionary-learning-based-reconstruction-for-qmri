function zoom_box = export_reconstruction_comparison(truth,reconstructions,mask,method_names,design,path,append)
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
    if nargin<7, append=false; end
    % Only fill missing options; explicitly supplied settings take precedence.
    if ~isfield(design,'ZoomBox'), design.ZoomBox=[]; end
    if ~isfield(design,'ZoomInsetFraction'), design.ZoomInsetFraction=.43; end
    if ~isfield(design,'ZoomColor'), design.ZoomColor=[1 .12 .28]; end
    if ~isfield(design,'ComparisonFigureSize'), design.ComparisonFigureSize=[1100 1320]; end
    [n,m,channels]=size(truth);
    % assert stops with an explanatory message if a required condition fails.
    % These checks prevent misleading comparisons of differently sized maps.
    assert(channels==3 && n>=2 && m>=2,'Truth must contain three 2-D parameter maps.');
    assert(isequal(size(mask),[n m]),'Mask and parameter maps must have the same size.');
    assert(numel(reconstructions)==numel(method_names) && ~isempty(reconstructions), ...
        'Provide one name per reconstruction.');
    for k=1:numel(reconstructions)
        assert(isequal(size(reconstructions{k}),size(truth)), ...
            'Each reconstruction must have the same size as truth.');
    end
    % Limit inset size and require valid RGB values and positive figure sizes.
    validateattributes(design.ZoomInsetFraction,{'numeric'}, ...
        {'real','finite','scalar','>=',.15,'<=',.65});
    validateattributes(design.ZoomColor,{'numeric'}, ...
        {'real','finite','vector','numel',3,'>=',0,'<=',1});
    validateattributes(design.ComparisonFigureSize,{'numeric'}, ...
        {'real','finite','vector','numel',2,'positive'});
    %% 2. Select one common zoom region for every map and error panel
    zoom_box=design.ZoomBox;
    if isempty(zoom_box)
        % Locate the body's bounding box; use the full image for an empty mask.
        [body_rows,body_cols]=find(mask);
        if isempty(body_rows), body_rows=[1;n]; body_cols=[1;m]; end
        row_extent=max(body_rows)-min(body_rows);
        col_extent=max(body_cols)-min(body_cols);
        % Choose a lower-left crop relative to that box, clipped to the image.
        % min/max also ensure that the crop contains at least two pixels per side.
        first_row=min(n-1,max(1,round(min(body_rows)+.56*row_extent)));
        last_row=min(n,max(first_row+1,round(min(body_rows)+.72*row_extent)));
        first_col=min(m-1,max(1,round(min(body_cols)+.22*col_extent)));
        last_col=min(m,max(first_col+1,round(min(body_cols)+.38*col_extent)));
        zoom_box=[first_row last_row first_col last_col];
    end
    % A manually supplied crop must use valid, inclusive integer pixel indices.
    validateattributes(zoom_box,{'numeric'}, ...
        {'real','finite','integer','vector','numel',4,'positive'});
    assert(zoom_box(1)<zoom_box(2) && zoom_box(2)<=n && ...
        zoom_box(3)<zoom_box(4) && zoom_box(4)<=m, ...
        'ZoomBox must select at least 2 rows/columns within the reduced map grid.');
    % Turn the four bounding indices into row/column vectors for direct indexing.
    rows=zoom_box(1):zoom_box(2); cols=zoom_box(3):zoom_box(4);

    %% 3. Arrange the column headings and create the comparison figure
    % Prepend the reference to the methods; its own error panel is therefore zero.
    shown=[{truth} reshape(reconstructions,1,[])];
    headings=[{'Ground truth'} reshape(method_names,1,[])];
    % Short publication headings for the demo's two dictionary methods.
    headings(strcmp(headings,'One-step DL'))={'DL (one step)'};
    headings(strcmp(headings,'Nested DL'))={'DL (nested)'};
    % Odd rows show maps; the following even rows show absolute errors.
    labels={'T_1 maps (ms)','T_1 absolute error (ms)', ...
        'T_2 maps (ms)','T_2 absolute error (ms)', ...
        '\rho maps','\rho absolute error'};
    % Hidden figures can still be exported for unattended demo runs.
    visible='off'; if design.ShowFigures, visible='on'; end
    fig=figure('Name','Map and error comparison with zoom insets','Color','w', ...
        'Visible',visible,'Units','pixels','Position',[40 40 design.ComparisonFigureSize], ...
        'DefaultAxesFontName',design.FontName,'DefaultAxesFontSize',design.FontSize, ...
        'DefaultTextFontName',design.FontName,'DefaultTextFontSize',design.FontSize);
    % exportgraphics crops to visible content; this white frame preserves
    % page margins around titles, rotated labels and the footer.
    annotation(fig,'rectangle',[.001 .001 .998 .998],'Color','w','LineWidth',.1);
    % Explicit positions keep ALL panels equally sized, including the ones
    % beside the shared colorbars. Compute in pixels to preserve the layout.
    %% 4. Compute a uniform six-row layout with space for labels and colorbars
    position=get(fig,'Position'); width=position(3); height=position(4);
    columns=numel(shown); gap_x=18; gap_y=16; group_gap=12;
    left=65; right=72; top=49; bottom=53;
    % Fit square panel positions into both available width and available height.
    % The extra group gaps separate the T1, T2 and density row pairs.
    side=min((width-left-right-(columns-1)*gap_x)/columns, ...
        (height-top-bottom-5*gap_y-2*group_gap)/6);
    assert(side>0,'ComparisonFigureSize is too small for the comparison grid.');
    % Center the complete grid horizontally within the reserved margins.
    left=left+(width-left-right-columns*side-(columns-1)*gap_x)/2;

    %% 5. Draw maps, absolute errors and identical zoom windows
    for row=1:6
        % Rows 1/2 use channel 1, rows 3/4 channel 2, and rows 5/6 channel 3.
        channel=ceil(row/2); is_error=mod(row,2)==0;
        % Select these once per row so all methods have identical color scaling.
        if is_error
            limits=design.ErrorLimits(channel,:); colors=design.ErrorColors;
        else
            limits=design.MapLimits(channel,:); colors=design.MapColors{channel};
        end
        % Figure coordinates start at the bottom; rows are filled from the top.
        bottom_y=height-top-row*side-(row-1)*gap_y-floor((row-1)/2)*group_gap;
        for column=1:columns
            % Convert the calculated pixel position into normalized figure units.
            pos=[left+(column-1)*(side+gap_x) bottom_y side side]./[width height width height];
            ax=axes('Parent',fig,'Units','normalized','Position',pos);
            % Work on a local display array; the stored reconstruction is unchanged.
            values=shown{column}(:,:,channel);
            % Error is in the same physical units as the corresponding parameter.
            if is_error, values=abs(values-truth(:,:,channel)); end
            if design.ZeroBackground, values(~mask)=0; end
            % Draw the full image, crop box, connector and magnified crop together.
            draw_panel(ax,values,rows,cols,limits,colors,design);
            % A single title at the top identifies each complete method column.
            if row==1
                title(ax,headings{column},'FontWeight','normal', ...
                    'FontSize',design.FontSize+3,'Interpreter','none');
            end
            % A single rotated label at the left identifies each complete row.
            if column==1
                text(ax,-.10,.5,labels{row},'Units','normalized','Rotation',90, ...
                    'HorizontalAlignment','center','VerticalAlignment','bottom', ...
                    'FontSize',design.FontSize+1,'Interpreter','tex','Clipping','off');
            end
            % Share the last panel's colorbar across its entire row.
            if column==columns
                cb=colorbar(ax);
                cb.Units='normalized';
                cb.Position=[pos(1)+pos(3)+9/width pos(2) 10/width pos(4)];
                cb.FontSize=design.FontSize-1;
                % Creating a colorbar can shrink axes; restore the common panel size.
                ax.Position=pos;
            end
        end
    end
    %% 6. Record the crop coordinates and write or append the PDF page
    % The footer records both the crop and whether the body mask was displayed.
    background='Background shown';
    if design.ZeroBackground, background='Background masked for display'; end
    caption=sprintf(['Zoom: rows %d-%d, columns %d-%d (reconstruction grid).  %s.\n' ...
        'Absolute errors in physical units; shared color limits within each row; values above the limits saturate.'], ...
        zoom_box,background);
    annotation(fig,'textbox',[.045 .006 .94 .03],'String',caption, ...
        'EdgeColor','none','FontName',design.FontName,'FontSize',design.FontSize-1, ...
        'Interpreter','none','VerticalAlignment','middle','HorizontalAlignment','center');
    % Finish graphics layout before exporting; APPEND preserves earlier PDF pages.
    drawnow;
    exportgraphics(fig,path,'ContentType','vector','BackgroundColor','white', ...
        'Resolution',design.Resolution,'Append',append);
    % Hidden figures have served their purpose once the PDF has been written.
    if ~design.ShowFigures, close(fig); end
end

function draw_panel(ax,values,rows,cols,limits,colors,design)
    %% Draw the full map and its inset in the same image coordinate system
    [n,m]=size(values);
    % imagesc maps physical values through the specified common row limits.
    imagesc(ax,values,limits); colormap(ax,colors);
    % Equal pixel aspect ratio preserves anatomy; hold keeps the image when adding lines.
    axis(ax,'image'); axis(ax,'off'); hold(ax,'on');
    if ~isempty(ax.Toolbar), ax.Toolbar.Visible='off'; end
    % Keep the crop's aspect ratio. Insets are images in the SAME axes as
    % their parent map, sharing its colormap and limits without resampling.
    crop=values(rows,cols); crop_size=size(crop);
    % Scale both crop dimensions equally; the longer side sets the requested size.
    magnification=design.ZoomInsetFraction*min(n,m)/max(crop_size);
    inset_width=crop_size(2)*magnification;
    inset_height=crop_size(1)*magnification;
    % Put the inset inside the upper-right image corner with a small margin.
    % Image row coordinates increase downward, so a small Y value is near the top.
    pad=.02*min(n,m);
    inset_left=m+.5-pad-inset_width; inset_top=.5+pad;
    % Pixel centers are integers, so subtract half a pixel to outline their edges.
    roi=[cols(1)-.5 rows(1)-.5 numel(cols) numel(rows)];
    % Connect the crop's upper-right corner to the inset's lower-left corner.
    line(ax,[roi(1)+roi(3) inset_left], ...
        [roi(2) inset_top+inset_height],'Color',design.ZoomColor,'LineWidth',.8);
    rectangle(ax,'Position',roi,'EdgeColor',design.ZoomColor,'LineWidth',.8);
    % XData/YData specify pixel CENTERS, whereas rectangles follow edges.
    image(ax,'CData',crop,'CDataMapping','scaled', ...
        'XData',inset_left+[.5 crop_size(2)-.5]*magnification, ...
        'YData',inset_top+[.5 crop_size(1)-.5]*magnification);
    % Outline the inset after drawing it, then add a subtle frame around the map.
    rectangle(ax,'Position',[inset_left inset_top inset_width inset_height], ...
        'EdgeColor',design.ZoomColor,'LineWidth',1);
    rectangle(ax,'Position',[.5 .5 m n],'EdgeColor',[.18 .18 .18],'LineWidth',.5);
    % Fix the original image extent so annotation objects cannot expand the axes.
    xlim(ax,[.5 m+.5]); ylim(ax,[.5 n+.5]);
end
