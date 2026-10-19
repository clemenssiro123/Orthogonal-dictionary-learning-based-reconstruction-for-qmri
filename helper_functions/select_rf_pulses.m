function [RF,indices]=select_rf_pulses(source_RF,L)
%SELECT_RF_PULSES Select L pulses across the full stored RF sequence.
%   [RF,indices]=select_rf_pulses(source_RF,L) returns row vectors containing
%   the selected pulses and their indices in source_RF. L must be an integer
%   from 1 to numel(source_RF). For L>=2, rounded, uniformly spaced indices
%   include both endpoints; L=1 selects the first pulse.
%
%   Source values (including complex phases) are preserved exactly. This is
%   selection in pulse-index space, without interpolation or filtering. It
%   does not select TR values or preserve the original acquisition duration;
%   the caller must supply L repetition times for the new pulse train.
%
%   Example: select_rf_pulses(1:9,3) returns [1 5 9].

    validateattributes(source_RF,{'numeric'},{'vector','nonempty','finite'}, ...
        mfilename,'source_RF',1);
    validateattributes(L,{'numeric'},{'scalar','real','finite','integer','positive'}, ...
        mfilename,'L',2);
    source_length=numel(source_RF);
    if L>source_length
        error('select_rf_pulses:TooManyPulses', ...
            'L=%g exceeds the %d stored RF pulses. Choose L between 1 and %d.', ...
            L,source_length,source_length);
    end
    if L==1
        indices=1;
    else
        indices=round(linspace(1,source_length,double(L)));
    end
    RF=reshape(source_RF(indices),1,[]);
end
