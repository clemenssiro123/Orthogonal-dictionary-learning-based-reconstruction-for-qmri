function [gx, gy] = dirichlet_gradient(u, h)
%DIRICHLET_GRADIENT Matrix-free differences including all zero boundary edges.
%   [GX,GY] = dirichlet_gradient(U,H) treats every N-by-M image pixel as
%   an interior unknown. Exterior nodes one spacing H beyond each of the
%   four sides are fixed to zero; outermost image pixels are NOT clamped.
%   GX is N-by-(M+1), GY is (N+1)-by-M; H defaults to 1.
%   Along a row, differences are [u(1), diff(u), -u(end)]/H.
%   No body mask is applied. Single/double and complex arrays are supported.
%
%   With [gx,gy] = dirichlet_gradient(u,h),
%   dirichlet_divergence(gx,gy,h) equals dirichlet_laplacian(u,h).
%   For alpha/2*(||gx||_F^2+||gy||_F^2), the Euclidean gradient is
%   -alpha*dirichlet_laplacian(u,h). Apply any quadrature factor separately.

    %% 1. Check image/spacing inputs and allocate the two edge grids
    if nargin < 2, h = 1; end
    % validateattributes stops with an error for unsupported image types or spacing.
    validateattributes(u, {'double','single'}, {'2d','nonempty','nonsparse'}, mfilename, 'u');
    validateattributes(h, {'numeric'}, {'real','finite','scalar','positive'}, mfilename, 'h');
    h = double(h);
    [n,m] = size(u);
    % Include one extra edge in each direction for the two exterior boundaries.
    % 'like',u preserves single/double precision and real/complex storage.
    gx = zeros(n,m+1,'like',u);
    gy = zeros(n+1,m,'like',u);

    %% 2. Horizontal differences: left boundary, interior edges, right boundary
    % First edge: image value minus the zero exterior value to its left.
    gx(:,1) = u(:,1);
    % Each interior edge joins neighboring image columns.
    gx(:,2:m) = u(:,2:m)-u(:,1:m-1);
    % Last edge: zero exterior value minus the last image column.
    gx(:,m+1) = -u(:,m);

    %% 3. Vertical differences: top boundary, interior edges, bottom boundary
    gy(1,:) = u(1,:);
    % Each interior edge joins neighboring image rows.
    gy(2:n,:) = u(2:n,:)-u(1:n-1,:);
    gy(n+1,:) = -u(n,:);

    %% 4. Convert value differences to derivatives per unit grid spacing
    gx = gx/h;
    gy = gy/h;
end
