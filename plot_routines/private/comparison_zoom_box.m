function [zoom_box,design] = comparison_zoom_box(mask,design) % Return one validated zoom region and fill missing inset display settings.
%COMPARISON_ZOOM_BOX Shared inset settings for the grid and compact PDF rows.
% MATLAB R2021b. Coordinates are inclusive indices on the reconstruction grid.
    if ~isfield(design,'ZoomBox'), design.ZoomBox=[]; end % Use automatic region selection when no ZoomBox setting was supplied.
    if ~isfield(design,'ZoomInsetFraction'), design.ZoomInsetFraction=.43; end % Default the longest inset side to 43 percent of the shorter map side.
    if ~isfield(design,'ZoomColor'), design.ZoomColor=[1 .12 .28]; end % Use a reddish border/connector color when no RGB setting was supplied.
    validateattributes(design.ZoomInsetFraction,{'numeric'}, ... % Require a numeric inset size setting.
        {'real','finite','scalar','>=',.15,'<=',.65}); % Accept only a finite scalar fraction between 0.15 and 0.65.
    validateattributes(design.ZoomColor,{'numeric'}, ... % Require a numeric color specification.
        {'real','finite','vector','numel',3,'>=',0,'<=',1}); % Accept exactly three finite RGB components, each between zero and one.
    [n,m]=size(mask); % Read the reconstruction grid dimensions from the mask.
    zoom_box=design.ZoomBox; % Start with the user's explicit box, or [] for automatic selection.
    if isempty(zoom_box) % Calculate a box only when the user did not specify one.
        [body_rows,body_cols]=find(mask); % Find all nonzero body-mask pixel coordinates.
        if isempty(body_rows), body_rows=[1;n]; body_cols=[1;m]; end % Use the full image bounds if the mask contains no foreground pixels.
        row_extent=max(body_rows)-min(body_rows); % Measure the body's bounding-box extent in the row direction.
        col_extent=max(body_cols)-min(body_cols); % Measure the body's bounding-box extent in the column direction.
        % Lower-left body region; at least two pixels in each direction.
        first_row=min(n-1,max(1,round(min(body_rows)+.56*row_extent))); % Place the first crop row 56 percent down the body, leaving room for a second row.
        last_row=min(n,max(first_row+1,round(min(body_rows)+.72*row_extent))); % End near 72 percent down the body, enforcing two rows and the image boundary.
        first_col=min(m-1,max(1,round(min(body_cols)+.22*col_extent))); % Place the first crop column 22 percent across the body, leaving room for a second column.
        last_col=min(m,max(first_col+1,round(min(body_cols)+.38*col_extent))); % End near 38 percent across the body, enforcing two columns and the image boundary.
        zoom_box=[first_row last_row first_col last_col]; % Store inclusive [first row, last row, first column, last column] indices.
    end
    validateattributes(zoom_box,{'numeric'}, ... % Check both automatically selected and user-supplied crop coordinates.
        {'real','finite','integer','vector','numel',4,'positive'}); % Require four positive finite integer indices.
    assert(zoom_box(1)<zoom_box(2) && zoom_box(2)<=n && ... % Require distinct ordered rows whose last index stays inside the image.
        zoom_box(3)<zoom_box(4) && zoom_box(4)<=m, ... % Require distinct ordered columns whose last index stays inside the image.
        'ZoomBox must select at least 2 rows/columns within the reduced map grid.'); % Explain the required minimum crop size if either bounds check fails.
end
