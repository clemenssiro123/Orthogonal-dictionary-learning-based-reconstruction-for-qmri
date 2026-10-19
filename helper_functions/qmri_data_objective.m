function [value, terms, gradient, Hdata, data_product] = qmri_data_objective(x, data, RF, TR, r, zeta, scales, h, reference_times)
%QMRI_DATA_OBJECTIVE Sampled Fourier least squares, without a spatial prior.
%   f(x)=h^2/2*||Q(u)-data||_2^2, u=[reference_times,0]+scales.*x.
%   Only measured rows enter the objective, gradient and exact GN product.
%   Hdata is the local h^2/r Bloch Gram approximation in scaled coordinates;
%   it is exact for full sampling (r=1). data_product applies the exact
%   sampled GN Hessian h^2*real(J'*J*V) to one or more real direction columns.
%   Both contain data curvature only; solvers add their own priors/damping.
%   x is N-by-M-by-3; data and B are N-by-M-by-L; P=N*M pixels per map.
%   RF/TR describe the pulse sequence; zeta(t):r:end are its measured rows.
%   Outputs 3/4/5 request the gradient, local Hessian and exact product handle.
%   Request only value/terms during backtracking to skip derivative work.

    %% 1. Convert scaled deviations to physical T1, T2 and density
    [n,m,~]=size(x); 
    P=n*m;
    % Reshape makes each reference/scale broadcast over its parameter map.
    u=reshape([reference_times(:).' 0],1,1,3)+x.*reshape(scales,1,1,3);
    times=u(:,:,1:2);
    assert(all(isfinite(times(:)) & times(:)>0), ...
        'Reference plus deviation must give strictly positive physical T1/T2.');
    rho=u(:,:,3);

    %% 2. Simulate magnetization and optional derivatives with respect to T1/T2
    if nargout>=3
        % B is the unit-density signal; d1=dB/dT1 and d2=dB/dT2.
        [B,d1,d2]=Blochmap_vectorized(RF,TR,u(:,:,1),u(:,:,2),0);
        % Accumulate one real derivative image for each unknown parameter.
        gradient_maps=zeros(size(x));
    else
        % An objective-only call needs magnetization but no derivatives.
        B=Blochmap_vectorized(RF,TR,u(:,:,1),u(:,:,2),0);
    end
    % Each pixel has a symmetric 3-by-3 local Gram matrix: six unique entries.
    if nargout>=4, gram=zeros(P,6); end

    %% 3. Accumulate measured-data error and its gradient over all frames
    value=0;
    for t=1:numel(TR)
        % Only these Cartesian Fourier rows were acquired in frame t.
        rows=zeta(t):r:n;
        % Density scales the Bloch signal; division by sqrt(P) makes FFT unitary.
        predicted=fft2(rho.*B(:,:,t))/sqrt(P);
        % Unmeasured entries stay zero, so they contribute neither error nor gradient.
        residual=complex(zeros(n,m));
        residual(rows,:)=predicted(rows,:)-data(rows,:,t);
        % Frobenius norm sums squared magnitudes; h^2 is the spatial quadrature weight.
        value=value+.5*h^2*norm(residual,'fro')^2;
        if nargout>=3
            % Adjoint of fft2(.)/sqrt(P): bring the residual back to image space.
            adjoint_residual=sqrt(P)*ifft2(residual);
            % Chain rule for scaled x: du_j/dx_j=scales(j).
            % Time sensitivities also contain rho; density sensitivity does not.
            q1=scales(1)*rho.*d1(:,:,t);
            q2=scales(2)*rho.*d2(:,:,t);
            q3=scales(3)*B(:,:,t);
            % Apply the conjugate sensitivities to the backprojected residual.
            % Parameters are real even though the measured signals are complex.
            gradient_maps(:,:,1)=gradient_maps(:,:,1)+h^2*real(conj(q1).*adjoint_residual);
            gradient_maps(:,:,2)=gradient_maps(:,:,2)+h^2*real(conj(q2).*adjoint_residual);
            gradient_maps(:,:,3)=gradient_maps(:,:,3)+h^2*real(conj(q3).*adjoint_residual);
            if nargout>=4
                % Flatten each sensitivity in MATLAB's column-major pixel order.
                q1=q1(:); q2=q2(:); q3=q3(:);
                % Store entries (11,12,13,22,23,33), summed over acquisition frames.
                gram=gram+[abs(q1).^2,real(conj(q1).*q2),real(conj(q1).*q3), ...
                    abs(q2).^2,real(conj(q2).*q3),abs(q3).^2];
            end
        end
    end

    %% 4. Return scalar terms and the stacked gradient
    terms.data=value;
    % Stack all T1 derivatives, then T2, then density, matching the solver's x(:).
    if nargout>=3, gradient=gradient_maps(:); end

    %% 5. Assemble the optional sparse, approximate data Hessian
    if nargout>=4
        blocks=cell(1,6);
        % spdiags puts one local Gram entry per pixel onto a diagonal block.
        % 1/r approximates uniform Fourier sampling; it is not the actual mask.
        for j=1:6, blocks{j}=spdiags(h^2/r*gram(:,j),0,P,P); end
        % Symmetric 3-by-3 blocks couple the three parameters at the same pixel.
        Hdata=[blocks{1} blocks{2} blocks{3}; ...
               blocks{2} blocks{4} blocks{5}; ...
               blocks{3} blocks{5} blocks{6}];
    end

    %% 6. Build the optional exact sampled Gauss-Newton product
    if nargout>=5
        % The handle holds sensitivities at this x for repeated QP/PCG products.
        % A new call at an accepted x creates a new handle; no persistent cache.
        sensitivities={scales(1)*rho.*d1,scales(2)*rho.*d2,scales(3)*B};
        % N-by-1-by-L mask broadcasts across the M image columns.
        sampling_mask=false(n,1,numel(TR));
        for t=1:numel(TR), sampling_mask(zeta(t):r:n,1,t)=true; end
        % Apply the sampled Fourier operator and its adjoint using FFTs.
        data_product=@(directions) qmri_gauss_newton_product(directions,sensitivities,sampling_mask,h);
    end
end
