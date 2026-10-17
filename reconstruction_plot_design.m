function design = reconstruction_plot_design(lower_bounds,upper_bounds,method_indices)
%RECONSTRUCTION_PLOT_DESIGN Plot defaults from demo_compare_all_reconstructions.
% MATLAB R2021b. Bounds set display scales only; no solver values are changed.
% METHOD_INDICES select colors/styles: 1 BLIP, 2 LM, 3 one-step, 4 nested.
% One pair of bounds also supports the standalone image-denoising demo.
    validateattributes(lower_bounds,{'numeric'},{'real','finite','vector'});
    validateattributes(upper_bounds,{'numeric'},{'real','finite','vector'});
    lower_bounds=lower_bounds(:); upper_bounds=upper_bounds(:);
    assert(numel(lower_bounds)==numel(upper_bounds) && ...
        all(upper_bounds>lower_bounds),'Provide increasing display bounds per map.');
    validateattributes(method_indices,{'numeric'}, ...
        {'vector','integer','>=',1,'<=',4});
    count=numel(upper_bounds);
    design.MapColors=repmat({turbo(256)},1,count);
    design.MapLimits=[min(zeros(count,1),lower_bounds) upper_bounds];
    design.ErrorColors=hot(256);
    design.ErrorLimits=[zeros(count,1) .2*(upper_bounds-lower_bounds)];
    colors=[.49 .18 .56; .47 .67 .19; .85 .33 .10; 0 .45 .74];
    styles={'-','--','-.','-'};
    design.MethodColors=colors(method_indices,:);
    design.LineStyles=styles(method_indices);
    design.FontName='Arial'; design.FontSize=10; design.LineWidth=1.4;
    design.ReportFontName='Consolas'; design.ReportFontSize=10;
    design.FigureSize=[1500 850]; design.Resolution=200;
    design.ShowFigures=true; design.ZeroBackground=true;
    design.ProfileRow=[]; design.ZoomBox=[];
    design.ZoomInsetFraction=.43; design.ZoomColor=[1 .12 .28];
    design.ComparisonFigureSize=[1100 1320];
    design.RowPanelSize=300; design.RowPanelGap=14;
    design.RowShowTitles=false; design.ShowRowFigures=false;
    design.MapNames={'T1 (ms)','T2 (ms)','Density'};
    if count==1
        design.MapNames={'Image'};
        design.ComparisonFigureSize=[1100 600];
    end
end
