function curve = draw_reconstruction_curve(x,y,method,design)
% Apply the same method color and line style to every report.
    % Keep a legend handle even when there are no completed updates.
    if isempty(x) && isempty(y), x=NaN; y=NaN; end
    curve=plot(x,y,'Color',design.MethodColors(method,:), ...
        'LineStyle',design.LineStyles{method},'LineWidth',design.LineWidth);
end
