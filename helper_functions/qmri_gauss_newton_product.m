function products = qmri_gauss_newton_product(directions, sensitivities, sampling_mask, h)
%QMRI_GAUSS_NEWTON_PRODUCT Exact sampled Fourier data GN, without a matrix.
%   The qMRI objective helpers supply the cached scaled Bloch
%   sensitivities Q={d(rho*B)/dx1,d(rho*B)/dx2,d(rho*B)/dx3}, each N-by-M-by-L.
%   DIRECTIONS is real, (3*N*M)-by-number_of_directions, in map-major order.
%   SAMPLING_MASK is logical N-by-1-by-L, broadcast across image columns.
%   Returns h^2*real(J'*J*DIRECTIONS), where J = S*unitaryFFT*Q.
%   The sqrt(N*M) factors of FFT and its adjoint cancel in J'*J.
%   Uses the actual mask, with no 1/r approximation. No Jacobian or dense
%   Hessian is assembled. Memory is O(N*M*L), not O((3*N*M)^2).
%   Several QP directions are handled one at a time to limit peak memory.

    %% 1. Prepare one output vector for each supplied direction
    [n,m,~] = size(sensitivities{1});
    products = zeros(size(directions));
    for column = 1:size(directions,2)
        % Undo map-major stacking: one image for each parameter's perturbation.
        v = reshape(directions(:,column),n,m,3);

        %% 2. Apply the Bloch signal Jacobian to this parameter perturbation
        % Linearized signal change = q_T1*v_T1 + q_T2*v_T2 + q_rho*v_rho.
        % Each N-by-M perturbation broadcasts over the L acquisition frames.
        signal = sensitivities{1}.*v(:,:,1)+sensitivities{2}.*v(:,:,2) ...
            +sensitivities{3}.*v(:,:,3);

        %% 3. Apply the sampled Fourier operator followed by its adjoint
        % fft2 acts on each frame; the mask zeros unmeasured Fourier rows.
        % ifft2 includes 1/(N*M), so the two unitary FFT scale factors cancel.
        backprojected = ifft2(sampling_mask.*fft2(signal));

        %% 4. Apply the adjoint Bloch Jacobian and stack the three result maps
        result = zeros(n,m,3);
        for channel = 1:3
            % Conjugation is the adjoint; sum(...,3) combines all frames.
            % real() reflects real unknowns, and h^2 is the data-fit grid weight.
            result(:,:,channel) = h^2*sum(real(conj(sensitivities{channel}).*backprojected),3);
        end
        % Match the original direction layout: T1 pixels, T2 pixels, density pixels.
        products(:,column) = result(:);
    end
end
