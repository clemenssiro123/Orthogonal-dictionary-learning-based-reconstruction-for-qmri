function lap = dirichlet_laplacian(u, h)
%DIRICHLET_LAPLACIAN Matrix-free five-point Laplacian, zero exterior values.
%   LAP = dirichlet_laplacian(U,H), H default 1, equals div(gradient(U)).
%   It uses the same boundaries as dirichlet_gradient, with negative sign:
%   <U,-LAP> = ||GX||_F^2 + ||GY||_F^2. Works for singleton dimensions too.
%   O(numel(U)) work/storage; does not allocate gradient arrays or matrices.
%   This is a LINEAR operator, also suitable for Hessian-vector products.

    %% 1. Check the image and positive grid spacing
    if nargin < 2, h = 1; end
    validateattributes(u, {'double','single'}, {'2d','nonempty','nonsparse'}, mfilename, 'u');
    validateattributes(h, {'numeric'}, {'real','finite','scalar','positive'}, mfilename, 'h');
    %% 2. Apply the five-point stencil: neighbors minus four times the center
    % Every pixel has central weight -4, including pixels at the image boundary.
    lap = -4*u;
    % Add the neighbor ABOVE wherever that neighbor is stored in the image.
    lap(2:end,:) = lap(2:end,:)+u(1:end-1,:);
    % Add the neighbor BELOW.
    lap(1:end-1,:) = lap(1:end-1,:)+u(2:end,:);
    % Add the neighbor to the LEFT.
    lap(:,2:end) = lap(:,2:end)+u(:,1:end-1);
    % Add the neighbor to the RIGHT. Missing exterior neighbors contribute zero.
    lap(:,1:end-1) = lap(:,1:end-1)+u(:,2:end);

    %% 3. Scale the second differences by the squared grid spacing
    lap = lap/double(h)^2;
end
