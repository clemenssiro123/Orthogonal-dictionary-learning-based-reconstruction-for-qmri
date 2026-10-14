function [bloch,diff_t1,diff_t2] = Blochmap_vectorized(RFpulses, TR ,T1, T2, df)
%
% This function calculates the MRF signal evolution for an IR-bSSFP-based readout
% for a single material (T1, T2, df) for the input RF pulse train and TR.
% and also calculate the derivatives of the Bloch mapping.
%
% This code is based on the single isochromat method popularized by Klaus
% Scheffler. This works for bSSFP as long as one can assume that both the
% flip angle and off-resonance frequency are essentially homogeneous inside
% a voxel.
%
% INPUTS:
% RFpulses Complex RF pulse train
% (1 x #timepoints)
%
% TR TR (1 x #timepoints)
%
% theta=[T1, T2], df T1, T2, and off-resonance frequency
% m1 the previous magnetization vector.
%
% Make sure these are all in the same
% time units!
%
% OUTPUTS:
%
% the single of the discrete Bloch map and its derivatives with respect to
% theta = [T1,T2].
%
%
% This code use part of the code by Dan Ma, Mark Griswold to implement
% the Bloch response.
% copy right: Guozhi Dong, Kostas Papafitsoros, Michael Hintermueller.
% April 2018 
%
% Vectorization by Clemens Sirotenko, February 2022
%
% Signal sampling is at TE = TR/2; the second half-TR relaxation follows.
% Off-resonance evolution is currently disabled, so df is ignored.

if nargin<5
    df = 0;
end

if nargin<4
    T2 = 50;
end

if nargin<3
    T1 = 1000;
end

if nargin<2
    TR = rand(1,1000)*4+10;
end

if nargin<1
    RFpulses = 10*pi./180*randn(1,1000);
end

[n1,n2] = size(T1);

T1_vec = T1(:);
T2_vec = T2(:);

indizes1 = T1_vec>0;
indizes2 = T2_vec>0;




N = length(TR);
D = length(T1_vec);  % We assume that T2 and df are the same size

A1_vec = zeros(D,1);
A2_vec = zeros(D,1);

A1_vec(indizes1) = 1./(T1_vec(indizes1));
A2_vec(indizes2) = 1./(T2_vec(indizes2));

rf = abs(RFpulses);
rph = angle(RFpulses);

dict = complex(zeros(D,N));
dict_DT1 = complex(zeros(D,N));
dict_DT2 = complex(zeros(D,N));


% Create and initialize the states matrix based on whether we invert or not

M1 = repmat([0 0 -1].',1,D);

dif_m_1 = repmat([0 0 0].',1,D);
dif_m_2 = repmat([0 0 0].',1,D);

%m1=[0 0 -1].';  % This assumes a perfect inversion pulse with no delay


for i=1:N
  
    % Define the linear dynamic operators for pulse i
    
    rx =    [   1.0 0.0 0.0;        % rotation matrix for pulse
                0.0 cos(rf(i)) sin(rf(i));
                0.0 -sin(rf(i)) cos(rf(i))];
            
    rdzp =  [   cos(rph(i)) sin(rph(i)) 0.0; % RF phase
                -sin(rph(i)) cos(rph(i)) 0.0;
                0.0 0.0 1.0];
    
    rdzm =  [   cos(-rph(i)) sin(-rph(i)) 0.0;
                -sin(-rph(i)) cos(-rph(i)) 0.0;
                0.0 0.0 1.0];
                 
    
    % relaxation terms

    %introduce these matrices to avoid division by zero and resulting nan's
    %most likely this is not the most clever way
    a1 = zeros(1,D);
    a2 = zeros(1,D);
    b1 = zeros(1,D);
    b2 = zeros(1,D);
    a1(indizes1) = exp(-TR(i)./(2*T1_vec(indizes1))');
    a2(indizes2) = exp(-TR(i)./(2*T2_vec(indizes2))');
    b2(indizes2) = TR(i)./(2*T2_vec(indizes2).^2)';
    b1(indizes1) = TR(i)./(2*T1_vec(indizes1).^2)';


    E12 = [a2; a2; a1];  %write the relaxation terms as a long matrix 
    

     
    A1 = [zeros(1,D);  zeros(1,D)  ; a1];
    B1 = [zeros(1,D);  zeros(1,D); b1];

    A2 = [a2;  a2; zeros(1,D)];
    B2 = [b2;  b2; zeros(1,D)];
    %diff T1
    
    
    dif_m_1 = (A1.*B1).*(rdzp*rx*rdzm*M1) + E12.*(rdzp*rx*rdzm*dif_m_1) - B1.*A1;
    dif_m_2 = (A2.*B2).*(rdzp*rx*rdzm*M1) + E12.*(rdzp*rx*rdzm*dif_m_2);
    

    M1 = rdzp*rx*rdzm*M1;   % do RF pulse for all states

    % relax 1st 1/2 of TR
    M1 = E12.*M1;
    M1(3,:) = M1(3,:) + (1 - a1);
   
    %beta = df.'.*TR(i)*2*pi; % beta is now a Dx1 vector associated with df

% Off-resonance for first 1/2 of TR
%     M1(1:2,:) = M1(1:2,:).*repmat([cos(beta./2); cos(beta/2)], 1, size(M1, 2)) + ...
%         M1(2:-1:1,:).*repmat([sin(beta./2); - sin(beta./2)], 1, size(M1, 2));
    
    
    %Sample assuming TE=TR/2     

    dict(:,i) = M1(1,:) + (1i).*M1(2,:);   
    dict_DT1(:,i) = dif_m_1(1,:) + (1i).*dif_m_1(2,:);
    dict_DT2(:,i) = dif_m_2(1,:) + (1i).*dif_m_2(2,:);
    
    %dict_DT1(indizes1,i) = dif_m_1(1,indizes1) + (1i).*dif_m_1(2,indizes1);
    %dict_DT2(indizes2,i) = dif_m_2(1,indizes2) + (1j).*dif_m_2(2,indizes2);

    % relaxation for the second half after TE

    dif_m_1 = (A1.*B1).*M1 + E12.*dif_m_1 - B1.*A1;
    dif_m_2 = (A2.*B2).*M1 + E12.*dif_m_2;

    M1 = E12.*M1;
    M1(3,:) = M1(3,:) + (1 - a1);

    % Off-resonance for 2nd 1/2 of TR (can be commented out if df = 0)
%     M1(1:2,:) = M1(1:2,:).*repmat([cos(beta./2);cos(beta/2)], 1, size(M1, 2)) + ...
%        M1(2:-1:1,:).*repmat([sin(beta./2);-sin(beta./2)], 1, size(M1, 2));
    

end 
bloch = reshape(dict,n1,n2,N);
diff_t1 = reshape(dict_DT1,n1,n2,N);
diff_t2 = reshape(dict_DT2,n1,n2,N);

end
