function [results,design] = comparison_plot_inputs(results,overrides)
% Load saved comparison data and merge optional display overrides.
% Reading a MAT file never starts a solver or modifies the saved results.
    if ischar(results) || (isstring(results) && isscalar(results))
        source=char(results);
        if isfolder(source), source=fullfile(source,'results.mat'); end
        results=load(source);
    end
    assert(isstruct(results) && isscalar(results),'Provide a saved results structure, MAT filename, or run folder.');
    if isfield(results,'design')
        design=results.design;
    else
        assert(isfield(results,'experiment'),'Results need design or experiment settings.');
        design=reconstruction_plot_design(results.experiment.LowerBounds, ...
            results.experiment.UpperBounds,1:numel(results.method_names));
    end
    if nargin>1 && ~isempty(overrides)
        assert(isstruct(overrides) && isscalar(overrides),'Design overrides must be a scalar structure.');
        for field=fieldnames(overrides).'
            design.(field{1})=overrides.(field{1});
        end
    end
    if ~isfield(design,'MapNames'), design.MapNames={'T1 (ms)','T2 (ms)','Density'}; end
end
