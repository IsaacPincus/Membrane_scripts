function fig_handle = total_lipids(xaxis_name, xscale, yaxis_name, yscale,...
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
colours_cb = ["#27E0D8", "#648FFF", "#1D0C6B", "#DC267F", "#FE6100", ...
    "#FFB000", "#016D24"];
plot(vs.param_1_vals, vs.S_B_vals./(1+vs.alpha_B_vals), '--', 'displayname', '$\frac{S_B}{1+\alpha_B}$')
r_phi_vals = vs.R_vals.*sin(vs.phi_vals);
no_bend_B_lip = (vs.d_vals.^2-r_phi_vals.^2*pi)./(1+vs.alpha_B_vals);
% plot(param_1_vals, no_bend_B_lip, 'k:', 'displayname', '$d^2-\pi r_\phi^2$')
plot(vs.param_1_vals, vs.S_A_vals./(1+vs.alpha_A_vals)+vs.S_B_vals./(1+vs.alpha_B_vals),...
    '-', 'displayname', '$\frac{S_A}{1+\alpha_A}+\frac{S_B}{1+\alpha_B}$ (total lipids)')
% yyaxis right
ylabel('Lipids')
% plot(init_vals, on*S_A_init, '-.', 'linewidth', 0.5, 'HandleVisibility','off');
plot(vs.param_1_vals, vs.S_A_vals./(1+vs.alpha_A_vals), '--', 'displayname', '$\frac{S_A}{1+\alpha_A}$')
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