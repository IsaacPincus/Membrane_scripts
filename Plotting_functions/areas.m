function fig_handle = areas(xaxis_name, xscale, yaxis_name, yscale,...
    vs, patch_lims, plot_curves)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

fig_handle = figure('Position',[400,100,800,600]);

hold on
axes1 = gca;
axes1.XScale = xscale;
axes1.YScale = yscale;
xlabel(xaxis_name)
ylabel(yaxis_name)
% plot(init_vals, on*S_B_init, '-.', 'linewidth', 0.5, 'HandleVisibility','off');
% plot(init_vals, on*(S_B_init+S_A_init), '-.', 'linewidth', 0.5, 'HandleVisibility','off');
plot(vs.param_1_vals, vs.S_B_vals, '--', 'displayname', '$S_B$')
r_phi_vals = vs.R_vals.*sin(vs.phi_vals);
plot(vs.param_1_vals, vs.d_vals.^2-r_phi_vals.^2*pi, 'k:', 'displayname', '$d^2-\pi r_\phi^2$')
plot(vs.param_1_vals, vs.S_A_vals+vs.S_B_vals, '-', 'displayname', '$S_\mathrm{total}$')
% yyaxis right
ylabel('Area')
% plot(init_vals, on*S_A_init, '-.', 'linewidth', 0.5, 'HandleVisibility','off');
plot(vs.param_1_vals, vs.S_A_vals, '--', 'displayname', '$S_A$')
if numel(patch_lims) ~= 1
    small_param = patch_lims(1);
    large_param = patch_lims(2);
    axes_ylims = get(gca, 'YLim');
    patch([small_param, small_param, large_param, large_param],...
          [axes_ylims(1), axes_ylims(2), axes_ylims(2), axes_ylims(1)], 'm', 'FaceAlpha', 0.05, ...
          'EdgeColor', 'None', 'HandleVisibility', 'off')
end
if numel(plot_curves) ~= 1
    plot_vert_lines(vs.param_1_vals(plot_curves), colours_cb)
end
legend


end