function draw_comparison_panel(ax,values,rows,cols,limits,colors,design) % Draw one map and its zoom inset; the calling exporter writes the PDF.
%DRAW_COMPARISON_PANEL Full map with a linked inset, without interpolation.
% Shared by the portrait report and the standalone comparison/initialization PDFs.
    [n,m]=size(values); % Read the full image height n and width m in pixels.
    imagesc(ax,values,limits); colormap(ax,colors); % Display values using the supplied limits, then assign this panel's colormap.
    axis(ax,'image'); axis(ax,'off'); hold(ax,'on'); % Keep square image pixels, hide axis decorations, and retain subsequent overlays.
    if ~isempty(ax.Toolbar), ax.Toolbar.Visible='off'; end % Hide the axes toolbar so interactive controls do not enter the exported plot.
    crop=values(rows,cols); crop_size=size(crop); % Extract the selected rows/columns and record the crop's height and width.
    magnification=design.ZoomInsetFraction*min(n,m)/max(crop_size); % Choose one scale factor for both crop dimensions to preserve its aspect ratio.
    inset_width=crop_size(2)*magnification; % Convert the crop's column count to its displayed inset width.
    inset_height=crop_size(1)*magnification; % Convert the crop's row count to its displayed inset height.
    pad=.02*min(n,m); % Leave a border equal to 2 percent of the shorter full-image dimension.
    inset_left=m+.5-pad-inset_width; inset_top=.5+pad; % Place the inset at the upper right; image row coordinates increase downward.
    % Pixel centers are integers; rectangles outline their half-pixel edges.
    roi=[cols(1)-.5 rows(1)-.5 numel(cols) numel(rows)]; % Describe the source region as [left top width height] along pixel edges.
    line(ax,[roi(1)+roi(3) inset_left], ... % Start a connector with x coordinates from the source box to the inset.
        [roi(2) inset_top+inset_height],'Color',design.ZoomColor,'LineWidth',.8); % Connect the source's upper-right corner to the inset's lower-left corner.
    rectangle(ax,'Position',roi,'EdgeColor',design.ZoomColor,'LineWidth',.8); % Outline the original region that is enlarged in the inset.
    image(ax,'CData',crop,'CDataMapping','scaled', ... % Draw the crop using the same axes color limits as the full image.
        'XData',inset_left+[.5 crop_size(2)-.5]*magnification, ... % Position the first/last inset pixel centers horizontally after magnification.
        'YData',inset_top+[.5 crop_size(1)-.5]*magnification); % Position the first/last inset pixel centers vertically after magnification.
    rectangle(ax,'Position',[inset_left inset_top inset_width inset_height], ... % Outline the inset using its outside edges rather than its pixel centers.
        'EdgeColor',design.ZoomColor,'LineWidth',1); % Apply the configured border color and a one-point line width.
    rectangle(ax,'Position',[.5 .5 m n],'EdgeColor',[.18 .18 .18],'LineWidth',.5); % Add a thin dark frame around the full image.
    xlim(ax,[.5 m+.5]); ylim(ax,[.5 n+.5]); % Lock both axes to the original image edges so overlays cannot enlarge them.
end
