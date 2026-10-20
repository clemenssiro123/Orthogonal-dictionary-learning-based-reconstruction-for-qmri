function sequence=extend_rf_sequence(original,target_length)
%EXTEND_RF_SEQUENCE Continue the supplied 205-frame acquisition to 1000.
%   sequence=extend_rf_sequence(original,1000) accepts the structure loaded
%   from data/sequence_parameter_middle_205.mat. It returns a structure to
%   save with save(filename,'-struct','sequence','-v7'); it writes no files.
%
%   The original RF train exactly repeats a 120-frame motif: 50 active
%   pulses, 10 zero-RF frames, the same 50 pulses at half amplitude, then
%   10 zero-RF frames. Continue this motif without interpolation or changing
%   its zero RF phase. RF magnitudes are radians (0 to about 64.34 degrees).
%   The first 205 RF values are preserved exactly.
%
%   The original TR train is not periodic over 120 frames. Repeat its full
%   205 values independently to supply one positive TR (ms) per RF entry
%   and preserve all original RF/TR pairs. The demos normally ignore this
%   stored TR train and use TR_ms*ones(1,L), with TR_ms=12 in both presets.
%   L and df retain their existing meanings. target_length is 205 to 1000.
%
%   This is a continuation of the project's simulation sequence, not an
%   optimized acquisition or an exact reproduction of a published sequence.
%   Alternating high/low sinusoidal flip-angle blocks and recovery gaps
%   follow the general MRF design in Ma et al., Nature 495 (2013), 187-192,
%   https://doi.org/10.1038/nature11971. The original project's RF values,
%   phases and TR values, rather than that paper's full protocol, are used.
%   MATLAB R2021b; no additional toolboxes.

    if nargin<2, target_length=1000; end
    validateattributes(original,{'struct'},{'scalar'},mfilename,'original',1);
    assert(all(isfield(original,{'RF','TR','L','df'})), ...
        'The source must contain RF, TR, L and df.');
    validateattributes(original.RF,{'double'}, ...
        {'vector','real','finite','nonnegative','numel',205});
    validateattributes(original.TR,{'double'}, ...
        {'vector','real','finite','positive','numel',205});
    assert(isequal(original.L,205),'Use the original 205-frame source file.');
    validateattributes(target_length,{'numeric'}, ...
        {'scalar','real','finite','integer','>=',205,'<=',1000});

    rf=reshape(original.RF,1,[]);
    tr=reshape(original.TR,1,[]);
    period=120;
    assert(isequal(rf,rf(mod(0:204,period)+1)), ...
        'The source RF does not have the expected 120-frame period.');
    assert(isequal(rf(61:110),rf(1:50)/2) && ...
        all(rf([51:60 111:120])==0), ...
        'The source RF must have alternating full/half blocks and recovery gaps.');

    sequence=original;
    sequence.L=double(target_length);
    sequence.RF=rf(mod(0:sequence.L-1,period)+1);
    sequence.TR=tr(mod(0:sequence.L-1,numel(tr))+1);
    sequence.RFDesign=struct( ...
        'Method','Continue original RF motif without interpolation', ...
        'Generator','helper_functions/extend_rf_sequence.m', ...
        'OriginalFile','sequence_parameter_middle_205.mat', ...
        'OriginalLength',205,'RFPeriod',period, ...
        'ActiveFramesPerBlock',50,'RecoveryFramesPerBlock',10, ...
        'AlternatingBlockScale',[1 0.5], ...
        'FlipAngleUnits','radians','PhaseConvention','Original zero RF phase', ...
        'TRUnits','ms','TRPolicy','Repeat original 205 TR values independently', ...
        'DefaultDemoTR_ms',12, ...
        'Reference','https://doi.org/10.1038/nature11971', ...
        'Scope','Simulation extension; not an optimized or published protocol');
end
