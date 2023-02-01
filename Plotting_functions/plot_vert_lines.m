function plot_vert_lines(x_vals, colours)
    for ii=1:length(x_vals)
        plot(ones(1,2)*x_vals(ii), [get(gca, 'YLim')],...
            ':o', 'linewidth', 0.5, 'HandleVisibility','off',...
            'Color',colours(ii), 'MarkerFaceColor',colours(ii));
    end
end