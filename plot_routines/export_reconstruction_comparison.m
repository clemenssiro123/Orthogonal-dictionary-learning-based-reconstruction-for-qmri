function zoom_box = export_reconstruction_comparison(truth,reconstructions,mask,method_names,design,path,append)
%EXPORT_RECONSTRUCTION_COMPARISON Write or append one map/error comparison.
% ZOOM_BOX = EXPORT_RECONSTRUCTION_COMPARISON(TRUTH,RECONSTRUCTIONS,MASK, ...
%     METHOD_NAMES,DESIGN,PATH,APPEND)
% TRUTH and each cell in RECONSTRUCTIONS contain N-by-M-by-3 physical maps
% [T1(ms),T2(ms),Density], or N-by-M single images. MASK is N-by-M;
% METHOD_NAMES follows the reconstruction cell order.
% DESIGN is a complete display structure. PATH names the destination PDF;
% its parent folder must exist. APPEND defaults to false (replace the file).
% True appends one page. ZOOM_BOX returns the inclusive crop actually used.
% ShowFigures=false closes the figure after export. Use PLOT_COMPARISON
% directly to obtain an editable figure without writing a PDF.
% MATLAB R2021b.
% See also PLOT_COMPARISON, EXPORT_RECONSTRUCTION_ROWS.
    if nargin<7, append=false; end
    results=struct('truth',truth,'reconstructions',{reconstructions}, ...
        'mask',mask,'method_names',{method_names},'design',design);
    [fig,zoom_box]=plot_comparison(results);
    export_reconstruction_pdf(fig,path,design,append);
end
