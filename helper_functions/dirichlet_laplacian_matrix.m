function LAP = dirichlet_laplacian_matrix(N, M, h)
%DIRICHLET_LAPLACIAN_MATRIX Sparse zero-Dirichlet Laplacian in U(:) ordering.
%   Vectorizing dirichlet_laplacian(U,H) gives LAP*U(:) for an N-by-M image.
%   H defaults to 1. Diagonal entries are -4/H^2, including at image edges.
%   Use this when a solver needs an explicit matrix (e.g. quadprog).
%   Only the most recent (N,M,H) matrix is cached, to reuse it in backtracking.
%   No gradients, vector-field operators, or Neumann operators are built.

    %% 1. Check dimensions/spacing and use floating-point arithmetic
    if nargin < 3, h = 1; end
    validateattributes(N, {'numeric'}, {'real','finite','scalar','integer','positive'}, mfilename, 'N');
    validateattributes(M, {'numeric'}, {'real','finite','scalar','integer','positive'}, mfilename, 'M');
    validateattributes(h, {'numeric'}, {'real','finite','scalar','positive'}, mfilename, 'h');
    N = double(N); M = double(M); h = double(h);

    %% 2. Reuse the last matrix when the grid geometry has not changed
    % persistent variables retain their values between calls to this function.
    persistent cached_size cached_h cached_lap
    if isequal(cached_size,[N M]) && isequal(cached_h,h)
        LAP = cached_lap;
        return;
    end

    %% 3. Assemble the 2-D Laplacian from two 1-D second-difference matrices
    % Ty acts along rows; Tx acts along columns.
    Ty = second_difference(N);
    Tx = second_difference(M);
    % MATLAB stacks columns in U(:): repeat Ty per column and use Tx to couple columns.
    % speye creates sparse identities; kron repeats/couples the corresponding blocks.
    LAP = (kron(speye(M),Ty)+kron(Tx,speye(N)))/h^2;
    % Save the geometry and matrix for the next call.
    cached_size = [N M]; cached_h = h; cached_lap = LAP;
end

%% Helper: one-dimensional second differences with zero exterior neighbors
function A = second_difference(n)
    % Put +1 on the first superdiagonal; n=1 correctly gives no neighbor entries.
    edge = sparse(1:n-1,2:n,ones(1,n-1),n,n);
    % -2 at each center and +1 at each existing neighbor; edge' adds the subdiagonal.
    A = spdiags(-2*ones(n,1),0,n,n)+edge+edge';
end
