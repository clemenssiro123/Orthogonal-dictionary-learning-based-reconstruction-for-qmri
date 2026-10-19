function [bloch,diff_t1,diff_t2] = Blochmap_vectorized(RFpulses,TR,T1,T2,df)
%BLOCHMAP_VECTORIZED IR-bSSFP signal and relaxation-time derivatives.
%   T1/T2 are equally sized maps; RFpulses/TR have one entry per frame.
%   RF magnitudes are flip angles in radians and RF phases set the rotation
%   axis. TR, T1 and T2 must use the same time unit (ms in this project).
%   Outputs have size [size(T1),numel(TR)] and assume unit proton density.
%
%   Magnetization starts at [0;0;-1] after a perfect inversion pulse.
%   Each frame applies RF rotation, half-TR relaxation, signal sampling at
%   TE=TR/2, and the remaining half-TR relaxation. Requesting derivatives
%   propagates dM/dT1 and dM/dT2 through these same operations.
%   DF is retained for call compatibility; off-resonance is not modeled.
%   Legacy convention: nonpositive relaxation times have zero decay factors
%   and zero decay-factor derivatives. Reconstruction solvers require T1/T2>0.
%
%   Based on the single-isochromat method of Klaus Scheffler and Bloch
%   response code by Dan Ma and Mark Griswold.
%   Copyright: Guozhi Dong, Kostas Papafitsoros, Michael Hintermueller, April 2018.
%   Vectorization: Clemens Sirotenko, February 2022.

    if nargin<4, T2=50; end
    if nargin<3, T1=1000; end
    if nargin<2, TR=rand(1,1000)*4+10; end
    if nargin<1, RFpulses=10*pi/180*randn(1,1000); end

    [n,m]=size(T1);
    T1=T1(:).'; T2=T2(:).';
    positive_t1=T1>0; positive_t2=T2>0;
    frames=numel(TR); pixels=numel(T1);
    flip=abs(RFpulses); phase=angle(RFpulses);
    signal=complex(zeros(pixels,frames));
    magnetization=repmat([0;0;-1],1,pixels);

    with_derivatives=nargout>=2;
    if with_derivatives
        signal_t1=complex(zeros(pixels,frames));
        signal_t2=complex(zeros(pixels,frames));
        derivative_t1=zeros(3,pixels);
        derivative_t2=zeros(3,pixels);
    end

    for t=1:frames
        % RF rotation: phase shift, rotation about x, then undo phase shift.
        rx=[1 0 0; 0 cos(flip(t)) sin(flip(t)); 0 -sin(flip(t)) cos(flip(t))];
        rz=[cos(phase(t)) sin(phase(t)) 0; -sin(phase(t)) cos(phase(t)) 0; 0 0 1];
        rotation=rz*rx*rz';

        a1=zeros(1,pixels); a2=zeros(1,pixels);
        a1(positive_t1)=exp(-TR(t)./(2*T1(positive_t1)));
        a2(positive_t2)=exp(-TR(t)./(2*T2(positive_t2)));
        decay=[a2;a2;a1];
        rotated=rotation*magnetization;

        if with_derivatives
            % d(exp(-TR/(2*T)))/dT = exp(-TR/(2*T))*TR/(2*T^2).
            b1=zeros(1,pixels); b2=zeros(1,pixels);
            b1(positive_t1)=TR(t)./(2*T1(positive_t1).^2);
            b2(positive_t2)=TR(t)./(2*T2(positive_t2).^2);
            decay_t1=[zeros(2,pixels);a1.*b1];
            decay_t2=[a2.*b2;a2.*b2;zeros(1,pixels)];
            % The -decay_t1 term differentiates longitudinal recovery (1-a1).
            derivative_t1=decay_t1.*rotated+decay.*(rotation*derivative_t1)-decay_t1;
            derivative_t2=decay_t2.*rotated+decay.*(rotation*derivative_t2);
        end

        % First half-TR relaxation, then sample the transverse signal.
        magnetization=decay.*rotated;
        magnetization(3,:)=magnetization(3,:)+(1-a1);
        signal(:,t)=magnetization(1,:)+1i*magnetization(2,:);
        if with_derivatives
            signal_t1(:,t)=derivative_t1(1,:)+1i*derivative_t1(2,:);
            signal_t2(:,t)=derivative_t2(1,:)+1i*derivative_t2(2,:);
            % Differentiate the remaining half-TR before updating the state.
            derivative_t1=decay_t1.*magnetization+decay.*derivative_t1-decay_t1;
            derivative_t2=decay_t2.*magnetization+decay.*derivative_t2;
        end
        magnetization=decay.*magnetization;
        magnetization(3,:)=magnetization(3,:)+(1-a1);
    end

    bloch=reshape(signal,n,m,frames);
    if with_derivatives
        diff_t1=reshape(signal_t1,n,m,frames);
        diff_t2=reshape(signal_t2,n,m,frames);
    end
end
