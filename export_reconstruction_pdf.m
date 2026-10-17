function export_reconstruction_pdf(fig,path,design,append)
%EXPORT_RECONSTRUCTION_PDF Write a report page using the comparison-demo style.
% MATLAB R2021b. Retry vector-renderer failures on a fresh hidden figure,
% then use an image page if needed. Commit the PDF only after a successful
% export, so a failed append cannot damage the preceding pages.
    if nargin<4, append=false; end
    assert(isgraphics(fig,'figure'),'The report figure must still exist.');
    if ~design.ShowFigures, cleanup=onCleanup(@() close_if_valid(fig)); end
    axes_list=findall(fig,'Type','axes');
    for ax=axes_list.'
        if ~isempty(ax.Toolbar), ax.Toolbar.Visible='off'; end
    end
    drawnow;
    folder=fileparts(path); if isempty(folder), folder=pwd; end
    working_path=[tempname(folder) '.pdf'];
    file_cleanup=onCleanup(@() delete_if_present(working_path));
    prepare_pdf(path,working_path,append);
    try
        write_page(fig,working_path,design,append,'vector');
    catch first_error
        if ~is_vector_error(first_error) || ~isgraphics(fig,'figure')
            rethrow(first_error);
        end
        % A fresh graphics tree avoids reusing the exporter's failed state.
        retry_fig=copyobj(fig,groot);
        retry_cleanup=onCleanup(@() close_if_valid(retry_fig));
        retry_fig.Visible='off'; drawnow;
        prepare_pdf(path,working_path,append);
        try
            write_page(retry_fig,working_path,design,append,'vector');
            warning('qmri:PDFVectorRetry', ...
                'Vector PDF export succeeded on a fresh figure: %s',path);
        catch retry_error
            if ~is_vector_error(retry_error), rethrow(retry_error); end
            prepare_pdf(path,working_path,append);
            write_page(retry_fig,working_path,design,append,'image');
            warning('qmri:PDFImageFallback', ...
                'Vector PDF export failed; this page was saved as an image at %g DPI: %s', ...
                design.Resolution,path);
        end
    end
    [ok,message]=movefile(working_path,path,'f');
    assert(ok,'Could not replace the report PDF: %s',message);
end

function write_page(fig,path,design,append,content_type)
    exportgraphics(fig,path,'ContentType',content_type,'BackgroundColor','white', ...
        'Resolution',design.Resolution,'Append',append);
end

function prepare_pdf(path,working_path,append)
    delete_if_present(working_path);
    if append && isfile(path)
        [ok,message]=copyfile(path,working_path);
        assert(ok,'Could not preserve the preceding PDF pages: %s',message);
    end
end

function tf=is_vector_error(err)
    % Do not conceal unrelated errors such as invalid options or file access.
    tf=contains(err.message,'OutputHelper');
end

function close_if_valid(fig)
    if isgraphics(fig,'figure'), close(fig); end
end

function delete_if_present(path)
    if isfile(path), delete(path); end
end
