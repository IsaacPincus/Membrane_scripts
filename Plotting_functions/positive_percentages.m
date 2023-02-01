function fig_handle = positive_percentages(xaxis_name, xscale, yaxis_name, yscale,...
    vs, patch_lims, plot_curves, curve_to_plot)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

fig_handle = figure('Position',[400,100,800,600]);

hold on
axes1 = gca;
axes1.XScale = xscale;
axes1.YScale = yscale;
xlabel(xaxis_name)
ylabel(yaxis_name)

lines = ["-", ":", ":", "--", ":", "--"];
colours = ['k', 'b', 'r', 'r', 'g', 'g'];
E_pos = vs.E_all(1,:)-vs.E_all(2,:);
E_pos_2 = sum(vs.E_all(3:end,:));
E_perc = vs.E_all./E_pos;
for ii=3:6
    plot(vs.param_1_vals, E_perc(ii,:), ...
        strcat(colours(ii),lines(ii)))
end
plot(vs.param_1_vals, E_perc(3,:)+E_perc(4,:), 'r-');
plot(vs.param_1_vals, E_perc(5,:)+E_perc(6,:), 'g-');
% axes_ylims = get(gca, 'YLim');
% patch([small_param, small_param, large_param, large_param],...
%       [axes_ylims(1), axes_ylims(2), axes_ylims(2), axes_ylims(1)],...
%       'm', 'FaceAlpha', 0.05, ...
%       'EdgeColor', 'None', 'HandleVisibility', 'off')
% plot_vert_lines(param_1_vals(plot_curves), colours_cb)
legend({'$E_\mathrm{stretch,A}$','$E_\mathrm{stretch,B}$',...
    '$E_\mathrm{bend,A}$','$E_\mathrm{bend,B}$',...
    '$E_\mathrm{stretch,Total}$','$E_\mathrm{bend,Total}$'},...
    'Box','off','location', 'best')
ylim([0,1])
end