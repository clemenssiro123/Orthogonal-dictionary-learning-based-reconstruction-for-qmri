function [sT1,sT2,inten]=BLIP(L,r,n1,n2,data,dict,zeta,T1set,T2set,mask)
%BLIP Legacy BLIP reconstruction for MR fingerprinting (Davies et al.).
% Each trial takes a data-fitting step in signal space and projects each
% pixel's time course onto a Bloch dictionary entry times a density value.
% The main reconstruction interface is now blip_reconstruction.m; this
% helper retains its original calling convention and iteration behavior.
%
% Algorithm reference:
% M. Davies, G. Puy, P. Vandergheynst and Y. Wiaux,
% "A Compressed Sensing Framework for Magnetic Resonance Fingerprinting,"
% SIAM Journal on Imaging Sciences, 7(4), 2623-2656, 2014.
% https://doi.org/10.1137/130947246
% Code provenance: legacy implementation; original software author,
% copyright holder and license have not been verified.
%
% INPUTS
% L: number of frames; [n1,n2]: image size; r: Cartesian sampling stride.
% data: [n1,n2,L] measured Fourier data; only sampled rows are used.
% dict: [number of entries,L] Bloch signals, one parameter pair per row.
% zeta(t): first sampled row in frame t; subsequent rows are r apart.
% T1set, T2set: lookup values, with T2 varying fastest in dictionary order.
% mask: [n1,n2] output mask; zero entries suppress the returned density.
% OUTPUTS: [n1,n2] T1, T2, and proton-density (inten) maps.

%% 1. Initialize the accepted signal, trial budget, and step size
% sig is the current accepted signal; newsig will hold a trial signal.
sig = zeros(n1,n2,L);
ite = 1;% First trial number.
ite_time = 20; % Maximum trials, including those rejected by the step check.
kappa = 0.99; % Safety factor in the allowable step-size estimate below.
mu = r; % Initial data-fitting step size, based on the sampling stride.
index=zeros(1,n1*n2); % One matched dictionary row for each image pixel.
% Store squared data-residual and sampled signal-step norms for each trial.
data_res_norm = zeros(ite_time,1);
sig_res_norm = data_res_norm;
newsig = zeros(size(sig));

%% 2. Take a data-fitting step and project onto the Bloch dictionary
while ite<=ite_time
    % Display the trial number and allocate [timepoints,pixels] for matching.
    disp(ite);
    exa = zeros(L,n1*n2);
    for t=1:L
        % Use a unitary Fourier transform: fft2 divided by sqrt(pixel count).
        tem = 1/sqrt(n1*n2)*fft2(sig(:,:,t));
        % Fill measured-minus-predicted residuals only on acquired rows.
        % Unsampled entries remain zero, so they do not enter the gradient.
        res = zeros(n1,n2);
        res(zeta(t):r:end,:) = data(zeta(t):r:end,:,t)-tem(zeta(t):r:end,:);
        data_res_norm(ite) = data_res_norm(ite) + norm(res(:))^2;    
        % The Fourier adjoint is sqrt(pixel count)*ifft2. Adding this
        % residual correction gives a negative-gradient step of length mu.
        newsig(:,:,t) = sig(:,:,t) + mu*sqrt(n1*n2)*ifft2(res);
        % Flatten this frame so each column of exa is one pixel's time course.
        exa(t,:) = reshape(newsig(:,:,t),[1,n1*n2]);
    end
    % Match normalized dictionary rows to each pixel using real correlation.
    [inner,index] = templatematch1(dict,exa);
    % Retrieve the original, unnormalized signals selected by that match.
    matchout = dict(index,:);
    % Convert correlation with a unit-norm row into the density multiplying
    % the original row. This legacy calculation does not clip negative values.
    inten = inner.'./sqrt(sum(real(matchout.*conj(matchout)),2));
    % Scale each selected time course by its pixel's estimated density.
    matchout = repmat(inten,1,L).*matchout;
    
    % Release temporary arrays; F_res is a retained legacy clear target.
    clear exa F_res tem res;

    % Fold [pixels,timepoints] back into an image sequence [n1,n2,L].
    newsig=reshape(matchout,[n1,n2,L]); 

    %% 3. Check the trial step against the sampling operator
    for t=1:L
        % Measure how strongly the proposed signal change affects the data.
        sig_res=1/sqrt(n1*n2)*fft2(newsig(:,:,t)-sig(:,:,t));
        sig_res=sig_res(zeta(t):r:end,:);
        sig_res_norm(ite)= sig_res_norm(ite)+norm(sig_res(:)).^2;
    end

    % Compare signal-space and measured-data step norms. kappa makes the
    % resulting upper bound on mu slightly more conservative.
    omega = kappa*norm(newsig(:)-sig(:))^2/sig_res_norm(ite);
     if mu>omega
        % Reject this trial: keep sig and halve mu for the next trial.
        mu=mu/2; 
     else
       % Accept the projected signal when the step-size check passes.
       sig=newsig;
         
     end
    % Rejected trials also use one of the 20 iterations in this legacy loop.
    ite=ite+1;
    clear sig_res matchout
end

%% 4. Convert the final trial's dictionary indices to T1 and T2 maps
% index and inten refer to the last trial, even if that trial was rejected;
% this is the existing behavior of this older implementation.
sT1 = zeros(n1,n2);
sT2 = sT1;
for i=1:length(index)
 % T2 varies fastest: quotient/remainder identify the T1 block and T2 entry.
 k=floor(index(i)/length(T2set));
 l=mod(index(i),length(T2set));
 if l==0
     % A zero remainder is the last T2 entry of T1 block k.
     sT1(i)=T1set(k);
     sT2(i)=T2set(end);
 else
     % Otherwise, MATLAB's one-based indexing places this entry in block k+1.
     sT1(i)=T1set(k+1);
     sT2(i)=T2set(l);
 end
end

%% 5. Reshape density and suppress output pixels outside the mask
inten=reshape(inten,n1,n2);
inten(mask==0)=0;
% Also suppress T1 and T2 wherever the resulting density is exactly zero.
sT1(inten==0)=0;
sT2(inten==0)=0;


