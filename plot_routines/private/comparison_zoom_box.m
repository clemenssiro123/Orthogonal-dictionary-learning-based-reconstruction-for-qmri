function [zoom_box,design] = comparison_zoom_box(mask,design)
%COMPARISON_ZOOM_BOX Shared inset settings for the grid and compact PDF rows.
% MATLAB R2021b. Coordinates are inclusive indices on the reconstruction grid.
    if ~isfield(design,'ZoomBox'), design.ZoomBox=[]; end
    if ~isfield(design,'ZoomInsetFraction'), design.ZoomInsetFraction=.43; end
    if ~isfield(design,'ZoomColor'), design.ZoomColor=[1 .12 .28]; end
    validateattributes(design.ZoomInsetFraction,{'numeric'}, ...
        {'real','finite','scalar','>=',.15,'<=',.65});
    validateattributes(design.ZoomColor,{'numeric'}, ...
        {'real','finite','vector','numel',3,'>=',0,'<=',1});
    [n,m]=size(mask);
    zoom_box=design.ZoomBox;
    if isempty(zoom_box)
        [body_rows,body_cols]=find(mask);
        if isempty(body_rows), body_rows=[1;n]; body_cols=[1;m]; end
        row_extent=max(body_rows)-min(body_rows);
        col_extent=max(body_cols)-min(body_cols);
        % Lower-left body region; at least two pixels in each direction.
        first_row=min(n-1,max(1,round(min(body_rows)+.56*row_extent)));
        last_row=min(n,max(first_row+1,round(min(body_rows)+.72*row_extent)));
        first_col=min(m-1,max(1,round(min(body_cols)+.22*col_extent)));
        last_col=min(m,max(first_col+1,round(min(body_cols)+.38*col_extent)));
        zoom_box=[first_row last_row first_col last_col];
    end
    validateattributes(zoom_box,{'numeric'}, ...
        {'real','finite','integer','vector','numel',4,'positive'});
    assert(zoom_box(1)<zoom_box(2) && zoom_box(2)<=n && ...
        zoom_box(3)<zoom_box(4) && zoom_box(4)<=m, ...
        'ZoomBox must select at least 2 rows/columns within the reduced map grid.');
end
