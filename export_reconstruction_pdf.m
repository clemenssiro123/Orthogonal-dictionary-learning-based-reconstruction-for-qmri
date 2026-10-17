function export_reconstruction_pdf(fig,path,design,append)
%EXPORT_RECONSTRUCTION_PDF Write a report page using the comparison-demo style.
% MATLAB R2021b. Hidden report figures are closed after export.
    if nargin<4, append=false; end
    if ~design.ShowFigures, cleanup=onCleanup(@() close(fig)); end
    axes_list=findall(fig,'Type','axes');
    for ax=axes_list.'
        if ~isempty(ax.Toolbar), ax.Toolbar.Visible='off'; end
    end
    drawnow;
    exportgraphics(fig,path,'ContentType','vector','BackgroundColor','white', ...
        'Resolution',design.Resolution,'Append',append);
end
