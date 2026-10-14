function [maxval,matchout] = templatematch1(dict,sig)
%TEMPLATEMATCH1 Match each pixel's signal to a row of an MRF dictionary.
% dict: [number of entries, timepoints], with one nonzero signal per row.
% sig:  [timepoints, number of pixels], with one measured signal per column.
% maxval:   [1, number of pixels], the best real-valued correlation.
% matchout: [1, number of pixels], the corresponding dictionary row index.
% Only dictionary rows are normalized; the signal amplitudes are preserved.
%
% The correlation matrix has size [number of entries, number of pixels].
% For large problems, the caller can split sig into batches of columns.
% Original implementation: Nov. 15, 2012, Dan Ma and Mark Griswold.

%% 1. Report missing inputs and normalize each dictionary signal
% This legacy check prints a message; it does not stop execution itself.
if nargin<2     
    disp('Error - 2 inputs are required!');      
end

% Dividing each row by its Euclidean norm makes the match independent of
% dictionary amplitude. Nonzero row norms are assumed by this helper.
for i=1:size(dict(:,:),1)
    dict(i,:)=dict(i,:)./norm(dict(i,:)); 
end

%% 2. Correlate every dictionary entry with every pixel signal
% These two matrix products give real(conj(dict)*sig). The score is the
% real part of the complex inner product, not its absolute magnitude.
innerproduct1 = real(dict(:,:))*real((sig(:,:)));
innerproduct2 = imag(dict(:,:))*imag((sig(:,:)));

% Add the real-real and imaginary-imaginary contributions for each pair.
innerproduct = innerproduct1 + innerproduct2;

%% 3. Return the largest score and its dictionary index for every pixel
% With multiple dictionary rows, max searches down each column. A tie is
% resolved by the first occurrence, following MATLAB's max convention.
[maxval,matchout] = max(innerproduct);
