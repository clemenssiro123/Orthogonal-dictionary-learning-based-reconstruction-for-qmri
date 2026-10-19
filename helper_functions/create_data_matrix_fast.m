function data_matrix = create_data_matrix_fast(data, patch_size, R_mat)
%CREATE_DATA_MATRIX_FAST Extract one vectorized image patch per column.
%   R_mat defines the patch locations and stride (see create_R_matrix).
%   Infer the number of patches from R_mat, so rectangular patches and
%   strides greater than one use the same extraction code.

    %% 1. Check the image, patch dimensions and extraction-operator size
    validateattributes(data, {'double'}, {'2d','nonempty'}, mfilename, 'data');
    validateattributes(patch_size, {'numeric'}, ...
        {'real','finite','vector','numel',2,'integer','positive'}, mfilename, 'patch_size');
    patch_area = prod(double(patch_size));
    % R needs one column per image pixel and patch_area rows per extracted patch.
    if size(R_mat,2) ~= numel(data) || size(R_mat,1) == 0 || ...
            mod(size(R_mat,1),patch_area) ~= 0
        error('create_data_matrix_fast:IncompatibleSize', ...
            'R_mat must have one column per image pixel and a whole number of patches.');
    end
    %% 2. Extract all patches and place one vectorized patch in each column
    % data(:) stacks image columns; R_mat copies entries in patch order.
    data_matrix = reshape(R_mat*data(:), patch_area, []);
end
