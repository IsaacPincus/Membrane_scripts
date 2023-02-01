function fig_handle = energy_per_lipid(xaxis_name, xscale, yaxis_name, yscale,...
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

lines = ["-", ":", ":", "--", ":", "--"];
colours = ['k', 'b', 'r', 'r', 'g', 'g'];
E_per_lipid(2,:) = -vs.E_all(2,:)./vs.S_A_vals.*(1+vs.alpha_A_vals);
E_per_lipid([3,5],:) = vs.E_all([3,5],:)./vs.S_A_vals.*(1+vs.alpha_A_vals);
E_per_lipid([4,6],:) = vs.E_all([4,6],:)./vs.S_B_vals.*(1+vs.alpha_B_vals);
for ii=3:6
    plot(vs.param_1_vals, E_per_lipid(ii,:), ...
        strcat(colours(ii),lines(ii)))
end
% plot(param_1_vals, E_per_lipid(3,:)+E_per_lipid(4,:), 'r-');
% plot(param_1_vals, E_per_lipid(5,:)+E_per_lipid(6,:), 'g-');
% legend({'$E_\mathrm{stretch,A}$','$E_\mathrm{stretch,B}$',...
%     '$E_\mathrm{bend,A}$','$E_\mathrm{bend,B}$',...
%     '$E_\mathrm{stretch,Total}$','$E_\mathrm{bend,Total}$'},...
%     'Box','off','location', 'best')
% legend({'$E_\mathrm{adhesion,A}$',...
%     '$E_\mathrm{stretch,A}$','$E_\mathrm{stretch,B}$',...
%     '$E_\mathrm{bend,A}$','$E_\mathrm{bend,B}$'}, 'Box','off',...
%     'location', 'best')
legend({'$E_\mathrm{stretch,A}$','$E_\mathrm{stretch,B}$',...
    '$E_\mathrm{bend,A}$','$E_\mathrm{bend,B}$'}, 'Box','off',...
    'location', 'best')

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