function div = dirichlet_divergence(gx, gy, h)
%DIRICHLET_DIVERGENCE Negative adjoint of dirichlet_gradient, without matrices.
%   GX must be N-by-(M+1), GY (N+1)-by-M; the result is N-by-M.
%   H defaults to 1. For real or complex arrays, <G*u,p> = -<u,div(p)>.
%   Consequently the gradient adjoint G'*p is -dirichlet_divergence(p).

    %% 1. Check edge-array types, spacing and compatible grid dimensions
    if nargin < 3, h = 1; end
    % Both components must represent edges of the SAME image grid.
    validateattributes(gx, {'double','single'}, {'2d','nonempty','nonsparse'}, mfilename, 'gx');
    validateattributes(gy, {'double','single'}, {'2d','nonempty','nonsparse'}, mfilename, 'gy');
    validateattributes(h, {'numeric'}, {'real','finite','scalar','positive'}, mfilename, 'h');
    if size(gy,1) ~= size(gx,1)+1 || size(gx,2) ~= size(gy,2)+1
        error('dirichlet_divergence:IncompatibleSize', ...
            'Expected gx of size N-by-(M+1) and gy of size (N+1)-by-M.');
    end
    %% 2. Sum right-minus-left and bottom-minus-top edge differences
    % The result has one value per image pixel; all exterior edges participate.
    % This sign convention makes divergence the NEGATIVE adjoint of gradient.
    div = (gx(:,2:end)-gx(:,1:end-1) + gy(2:end,:)-gy(1:end-1,:))/double(h);
end
