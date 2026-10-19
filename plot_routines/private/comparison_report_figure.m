function fig=comparison_report_figure(name,design) % Create a figure for the existing multi-page reports; export is handled separately.
    % Apply the same geometry/fonts to each report page; off still allows PDF export.
    visible='off'; if design.ShowFigures, visible='on'; end % Hide report windows unless ShowFigures requests an on-screen display.
    fig=figure('Name',name,'Color','w','Visible',visible, ... % Create a named white figure with the chosen visibility.
        'Position',[40 40 design.FigureSize],'DefaultAxesFontName',design.FontName, ... % Apply the configured screen position/size and axes font family.
        'DefaultAxesFontSize',design.FontSize,'DefaultTextFontName',design.FontName, ... % Set the axes font size and text font family for later plot elements.
        'DefaultTextFontSize',design.FontSize); % Set the default text font size for the new report figure.
end
