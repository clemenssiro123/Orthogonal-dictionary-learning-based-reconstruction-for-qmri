function paths = export_comparison_reports(results,output_folder,design)
%EXPORT_COMPARISON_REPORTS Recreate all four reports and eight compact PDFs.
% PATHS = EXPORT_COMPARISON_REPORTS(RESULTS,OUTPUT_FOLDER,DESIGN)
% RESULTS is the four-method results structure, MAT filename, or run folder.
% OUTPUT_FOLDER is created if needed. Existing PDFs with the same names in
% that folder are replaced; results.mat is not written by this function.
% Optional DESIGN overrides selected saved display fields. ShowFigures=false
% exports without visible windows and closes the generated hidden figures.
% ShowRowFigures additionally controls visibility of the compact rows.
%
% PATHS is a cell array of twelve filenames: maps, errors, convergence,
% profiles/summary, then eight compact map/error PDFs. The profile report
% includes an extra enlarged-map page only when ZoomBox is explicitly set.
% Page order and display settings match demo_compare_all_reconstructions.
% This function never reruns simulation or reconstruction. MATLAB R2021b.
%
% Example: export_comparison_reports(s,fullfile(run_folder,'replots'), ...
%              struct('ShowFigures',false));
% See also PLOT_COMPARISON, PLOT_CONVERGENCE, EXPORT_RECONSTRUCTION_PDF.
    if nargin<3, design=[]; end
    [results,design]=comparison_plot_inputs(results,design);
    if ~isfolder(output_folder), mkdir(output_folder); end
    maps_pdf=fullfile(output_folder,'maps.pdf');
    errors_pdf=fullfile(output_folder,'errors.pdf');
    convergence_pdf=fullfile(output_folder,'convergence.pdf');
    profiles_pdf=fullfile(output_folder,'profiles_and_summary.pdf');
    % Start each report with append=false, then add its remaining pages in order.
    export_reconstruction_pdf(plot_parameter_maps(results,design),maps_pdf,design);
    export_reconstruction_pdf(plot_comparison(results,design),maps_pdf,design,true);
    compact_paths=export_reconstruction_rows(results.truth,results.reconstructions, ...
        results.u0,results.mask,results.method_names,design,output_folder);
    export_reconstruction_pdf(plot_parameter_errors(results,design),errors_pdf,design);
    export_reconstruction_pdf(plot_convergence(results,'fit',design),convergence_pdf,design);
    export_reconstruction_pdf(plot_convergence(results,'errors',design),convergence_pdf,design,true);
    export_reconstruction_pdf(plot_convergence(results,'dictionary',design),convergence_pdf,design,true);
    export_reconstruction_pdf(plot_parameter_profiles(results,design),profiles_pdf,design);
    if ~isempty(design.ZoomBox)
        % Match the demo: an explicitly selected crop gets a separate page.
        export_reconstruction_pdf(plot_parameter_zoom(results,design),profiles_pdf,design,true);
    end
    export_reconstruction_pdf(plot_experiment_summary(results,design),profiles_pdf,design,true);
    paths=[{maps_pdf,errors_pdf,convergence_pdf,profiles_pdf} compact_paths];
end
