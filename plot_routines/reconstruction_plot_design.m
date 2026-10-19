function design = reconstruction_plot_design(lower_bounds,upper_bounds,method_indices)
%RECONSTRUCTION_PLOT_DESIGN Complete display defaults for standalone demos.
% DESIGN = RECONSTRUCTION_PLOT_DESIGN(LOWER_BOUNDS,UPPER_BOUNDS,METHOD_INDICES)
% Supply three bounds per side in physical [T1,T2,Density] units, or one
% bound pair for a single-image denoising plot. Bounds set display ranges;
% they do not alter solver constraints or stored reconstruction values.
% METHOD_INDICES selects method colors/styles in the requested order:
% 1 BLIP, 2 LM, 3 one-step DL, 4 nested DL.
% DESIGN includes map/error limits and colors, fonts, figure sizes, zoom,
% profile row and export visibility. Reuse a saved run's design to preserve
% any custom styling from that run instead of constructing fresh defaults.
% MATLAB R2021b.
% See also RECONSTRUCTION_REPORT_FIGURE, PLOT_COMPARISON.
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
