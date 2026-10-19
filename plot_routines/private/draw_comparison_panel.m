function draw_comparison_panel(ax,values,rows,cols,limits,colors,design)
%DRAW_COMPARISON_PANEL Full map with a linked inset, without interpolation.
% Shared by the portrait report and the standalone comparison/initialization PDFs.
    [n,m]=size(values);
    imagesc(ax,values,limits); colormap(ax,colors);
    axis(ax,'image'); axis(ax,'off'); hold(ax,'on');
    if ~isempty(ax.Toolbar), ax.Toolbar.Visible='off'; end
    crop=values(rows,cols); crop_size=size(crop);
    % Scale both crop dimensions equally and place the inset at the upper right.
    magnification=design.ZoomInsetFraction*min(n,m)/max(crop_size);
    inset_width=crop_size(2)*magnification;
    inset_height=crop_size(1)*magnification;
    pad=.02*min(n,m);
    inset_left=m+.5-pad-inset_width; inset_top=.5+pad;
    % Pixel centers are integers; rectangles outline their half-pixel edges.
    roi=[cols(1)-.5 rows(1)-.5 numel(cols) numel(rows)];
    line(ax,[roi(1)+roi(3) inset_left], ...
        [roi(2) inset_top+inset_height],'Color',design.ZoomColor,'LineWidth',.8);
    rectangle(ax,'Position',roi,'EdgeColor',design.ZoomColor,'LineWidth',.8);
    image(ax,'CData',crop,'CDataMapping','scaled', ...
        'XData',inset_left+[.5 crop_size(2)-.5]*magnification, ...
        'YData',inset_top+[.5 crop_size(1)-.5]*magnification);
    rectangle(ax,'Position',[inset_left inset_top inset_width inset_height], ...
        'EdgeColor',design.ZoomColor,'LineWidth',1);
    rectangle(ax,'Position',[.5 .5 m n],'EdgeColor',[.18 .18 .18],'LineWidth',.5);
    % Overlays must not change the displayed extent of the full image.
    xlim(ax,[.5 m+.5]); ylim(ax,[.5 n+.5]);
end
