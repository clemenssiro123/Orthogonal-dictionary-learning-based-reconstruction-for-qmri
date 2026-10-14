function [R, number_of_patches, overlap_count] = create_R_matrix(image_size, patch_size, r)
%CREATE_R_MATRIX Sparse extraction operator for all valid patches on a grid.
%   [R,M,D] = create_R_matrix([N1 N2],[P1 P2],RSTEP) uses stride RSTEP
%   in both image directions. Sizes and stride must be positive integers;
%   the patch must fit inside the image. Rectangular images/patches are valid.
%
%   For an N1-by-N2 image U, reshape(R*U(:),P1*P2,M) contains one patch
%   per column. Pixels within each patch follow MATLAB column-major order.
%   Patch origins visit rows first, then columns, starting at (1,1).
%   Only fully contained patches are used: no padding or periodic wrapping.
%   Origins are 1:RSTEP:(Ni-Pi+1); the last origin need not touch the edge.
%
%   R has (P1*P2*M)-by-(N1*N2) entries and one nonzero (1) per row.
%   Optional D is the N1*N2-by-1 vector of pixel overlap counts, satisfying
%   R'*R = spdiags(D,0,N1*N2,N1*N2). Compute/cache it once per geometry.
%   R' sums patch values back into the image; divide by D to average them.
%   For stride > 1 some pixels may be uncovered (D=0), so do not divide
%   those entries by zero. Stride 1 guarantees coverage for valid sizes.
%
%   Existing two-output calls remain supported. MATLAB R2021b compatible.

    %% 1. Check the image size, patch size and stride
    % validateattributes rejects malformed inputs before any indexing is attempted.
    validateattributes(image_size, {'numeric'}, ...
        {'real','finite','vector','numel',2,'integer','positive'}, mfilename, 'image_size');
    validateattributes(patch_size, {'numeric'}, ...
        {'real','finite','vector','numel',2,'integer','positive'}, mfilename, 'patch_size');
    validateattributes(r, {'numeric'}, ...
        {'real','finite','scalar','integer','positive'}, mfilename, 'r');
    % Normalize row/column vectors and integer-valued numeric input classes.
    image_size = double(image_size(:));
    patch_size = double(patch_size(:));
    r = double(r);
    if any(patch_size > image_size)
        error('create_R_matrix:PatchTooLarge', 'The patch must fit inside the image.');
    end

    % Image dimensions n1/n2 and patch dimensions p1/p2 are rows/columns.
    n1 = image_size(1);
    n2 = image_size(2);
    p1 = patch_size(1);
    p2 = patch_size(2);

    %% 2. Enumerate all valid top-left patch origins
    % ndgrid combines every allowed row with every allowed column, without padding.
    [start_rows, start_columns] = ndgrid(1:r:(n1-p1+1), 1:r:(n2-p2+1));
    % A pixel (row,col) has linear image index row+n1*(col-1).
    origins = start_rows(:) + n1*(start_columns(:)-1);
    number_of_patches = numel(origins);

    %% 3. List the image pixels belonging to every patch
    % Relative pixel indices within one patch, in the same order as U(:).
    offsets = (0:p1-1).' + n1*(0:p2-1);
    % Implicit expansion adds every within-patch offset to every patch origin.
    % Each column now lists the image indices for one patch.
    pixel_indices = offsets(:) + origins.';
    entries = numel(pixel_indices);

    %% 4. Make a sparse selector with one unit entry per extracted patch pixel
    % Each row copies one image pixel. Repeated columns represent overlapping patches.
    R = sparse((1:entries).', pixel_indices(:), ones(entries,1), entries, n1*n2);

    %% 5. Optionally count how often each image pixel is extracted
    if nargout >= 3
        % Column sums give diag(R'*R), without multiplying the two sparse matrices.
        overlap_count = full(sum(R,1)).';
    end
end
