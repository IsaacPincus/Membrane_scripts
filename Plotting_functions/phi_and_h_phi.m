function fig_handle = phi_and_h_phi(xaxis_name, xscale, yaxis_name, yscale,...
    vs, patch_lims, plot_curves,fig_handle_inp)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

colours_cb = ["#27E0D8", "#648FFF", "#1D0C6B", "#DC267F", "#FE6100", ...
    "#FFB000", "#016D24"];

if exist('fig_handle_inp', 'var')
    figure(fig_handle_inp)
    fig_handle = fig_handle_inp;
else
    fig_handle = figure('Position',[400,100,800,600]);
end
hold on
axes1 = gca;
axes1.XScale = xscale;
axes1.YScale = yscale;
xlabel(xaxis_name)
ylabel(yaxis_name)

% plot(vs.param_1_vals, rad2deg(vs.phi_vals), 'displayname', '$\phi$')
plot(vs.param_1_vals, rad2deg(vs.phi_vals),...
    'displayname', sprintf('$R=%0.2g$ $\\mu$m, $\\zeta=%0.2g$, $\\alpha_i=%0.2g$', vs.R, vs.zeta, vs.alpha_i))

% plot(init_vals, on*rad2deg(phi_init), '-.', 'linewidth', 0.5, 'HandleVisibility','off')
% yyaxis right
% ylabel('$h_\phi$')
% plot(param_1_vals, h_phi_vals, 'displayname', '$h_\phi$')

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