function show_reconstruction_map(ax,values,mask,limits,colors,design)
% Mask a display copy; stored physical values are never altered.
    if design.ZeroBackground, values(~mask)=0; end
    imagesc(ax,values,limits); axis(ax,'image'); axis(ax,'off'); colormap(ax,colors);
end
