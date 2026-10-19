function prox = prox_c(gamma,c,mode)
%PROX_C Minimize a coefficient penalty plus squared distance to the input.
% For a nonnegative scalar gamma, each output entry independently solves
% min_z 0.5*|z-c|^2 + gamma*penalty(z). The array shape is preserved.
% 'soft' uses penalty(z)=|z|; the other branch uses penalty(z)=(z~=0).
% This helper is used by dictionary learning to make small coefficients zero.

%% 1. Soft thresholding: shrink every magnitude by gamma, stopping at zero
if strcmp(mode,'soft')
    % abs measures size; sign restores the original sign (or complex phase).
    % Magnitudes <= gamma vanish; larger magnitudes are shortened by gamma.
    prox = sign(c).*max(abs(c) - gamma,0);
else
    %% 2. Hard thresholding: retain large coefficients without shrinking them
    % Keeping c costs gamma; replacing it by zero costs |c|^2/2.
    prox = c;
    % The two choices have equal cost at |c|=sqrt(2*gamma).
    alpha = sqrt(2*gamma);
    % At the tie this implementation selects zero. Every non-'soft' mode
    % enters this branch; the intended mode name is 'hard'.
    prox(abs(c)<=alpha) = 0;
end
end
