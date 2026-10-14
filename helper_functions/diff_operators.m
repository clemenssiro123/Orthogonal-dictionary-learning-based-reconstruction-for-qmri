function [DX,DY,div_X,div_Y,LAP] = diff_operators(N,M,h)
%DIFF_OPERATORS Consistent sparse scalar operators with zero Dirichlet exterior.
%   [DX,DY,DIV_X,DIV_Y,LAP] = diff_operators(N,M,H), H default 1.
%   All N-by-M image pixels are interior unknowns, with zero exterior nodes
%   one spacing H beyond all four sides. No body support mask is imposed.
%
%   GX = reshape(DX*u(:),N,M+1), GY = reshape(DY*u(:),N+1,M).
%   Gradient edges include BOTH exterior boundaries in each direction.
%   DIV_X = -DX', DIV_Y = -DY'; LAP = DIV_X*DX + DIV_Y*DY.
%   Therefore -u'*LAP*u = norm(DX*u)^2 + norm(DY*u)^2.
%
%   API CHANGE: LAP is output 5. The unused legacy GRADDIV, E and Neumann
%   lapn outputs were removed. All reconstruction callers use the dedicated
%   dirichlet_laplacian_matrix or matrix-free dirichlet_laplacian instead.
%   For matrix-free gradient/divergence, see dirichlet_gradient/divergence.
%   This routine builds matrices for inspection or algorithms that need them;
%   do not call it merely to apply a Laplacian in a repeated inner loop.

    %% 1. Validate grid dimensions and spacing
    if nargin < 3, h = 1; end
    validateattributes(N, {'numeric'}, {'real','finite','scalar','integer','positive'}, mfilename, 'N');
    validateattributes(M, {'numeric'}, {'real','finite','scalar','integer','positive'}, mfilename, 'M');
    validateattributes(h, {'numeric'}, {'real','finite','scalar','positive'}, mfilename, 'h');
    N = double(N); M = double(M); h = double(h);

    %% 2. Build directional gradients in MATLAB's column-major image ordering
    % A line with n pixels has n+1 edges, including both zero-exterior boundaries.
    Gy = edge_difference(N)/h;
    Gx = edge_difference(M)/h;
    % DX couples columns; the identity preserves row positions inside each column.
    DX = kron(Gx,speye(N));
    % DY applies the row-direction difference independently in every image column.
    DY = kron(speye(M),Gy);

    %% 3. Define divergence as the negative adjoint of the gradient
    div_X = -DX';
    div_Y = -DY';

    %% 4. Build the Laplacian only when the caller requests the fifth output
    if nargout >= 5
        % Equivalent to div_X*DX + div_Y*DY, with the same exterior conditions.
        LAP = dirichlet_laplacian_matrix(N,M,h);
    end
end

%% Helper: sparse first differences along one line of pixels
function G = edge_difference(n)
% For one line: G*u = [u(1); diff(u); -u(end)].
    % sparse(row,column,value,rows,cols) places +u_i in edge i and -u_i in edge i+1.
    G = sparse([1:n 2:n+1],[1:n 1:n],[ones(1,n) -ones(1,n)],n+1,n);
end
