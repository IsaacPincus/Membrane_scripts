function fig_handle = delta_E(xaxis_name, xscale, yaxis_name, yscale,...
    vs, patch_lims, plot_curves, anno_string)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

param_1_vals = vs.param_1_vals;
E_all = vs.E_all.*vs.sigma_vals;
init_vals = vs.alpha_i^2*vs.kD*vs.d_vals.^2/2.*vs.sigma_vals;
% E_all = vs.E_all;
% init_vals = vs.alpha_i^2*vs.kD*vs.d_vals.^2/2;

fig_handle = figure('Position',[400,100,800,600]);

hold on
axes1 = gca;
axes1.XScale = xscale;
axes1.YScale = yscale;
xlabel(xaxis_name)
ylabel(yaxis_name)
lines = ["-", ":", ":", "--", ":", "--"];
colours = ['k', 'b', 'r', 'r', 'g', 'g'];
colours_cb = ["#27E0D8", "#648FFF", "#1D0C6B", "#DC267F", "#FE6100", ...
    "#FFB000", "#016D24"];
E_all(1,:) = E_all(1,:) - init_vals;
E_all(4,:) = E_all(4,:) - init_vals;
for ii=1:6
%     plot(param_1_vals, E_all(ii,:)./(vs.kD.*vs.d_vals.^2), ...
%         strcat(colours(ii),lines(ii)))
    plot(param_1_vals, E_all(ii,:), ...
        strcat(colours(ii),lines(ii)))
end
min_E = min(E_all, [], 'all')*0.9;
max_E = max(E_all, [], 'all')*1.1;
if numel(patch_lims) ~= 1
    small_param = patch_lims(1);
    large_param = patch_lims(2);
    patch([small_param, small_param, large_param, large_param],...
          [min_E, max_E, max_E, min_E], 'm', 'FaceAlpha', 0.05, ...
          'EdgeColor', 'None', 'HandleVisibility', 'off')
end
if numel(plot_curves) ~= 1
    plot_vert_lines(param_1_vals(plot_curves), colours_cb)
end
legend({'$E_\mathrm{total}$','$E_\mathrm{adhesion}$',...
    '$E_\mathrm{stretch,A}$','$E_\mathrm{stretch,B}$',...
    '$E_\mathrm{bend,A}$','$E_\mathrm{bend,B}$'}, 'Box','off',...
    'location', 'best')
annotation('textbox', [0.3,0.6,0.4,0.2], 'String',...
    anno_string,...
    'interpreter', 'latex','FontSize',16,...
    'FitBoxToText','on','LineStyle','none', 'Color','b')
xlim([param_1_vals(1), param_1_vals(end)])
% ylim([min_E, max_E])
end