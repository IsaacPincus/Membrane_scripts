function fig_handle = delta_E_two_axes(xaxis_name, xscale, yaxis_name, yscale,...
    vs, D, patch_lims, plot_curves, anno_string)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

E_all = vs.E_all;
sigma_vals = vs.sigma_vals;
N_vals = sigma_vals*D^2./(pi*vs.R_vals);

fig_handle = figure('Position',[400,100,800,600]);
% t = tiledlayout(1,1);

hold on
% axes1 = axes(t);
axes1 = gca;
hold(axes1, 'on');
axes1.XScale = xscale;
axes1.YScale = yscale;
xlabel(xaxis_name)
ylabel(yaxis_name)
lines = ["-", ":", ":", "--", ":", "--"];
colours = ['k', 'b', 'r', 'r', 'g', 'g'];
colours_cb = ["#27E0D8", "#648FFF", "#1D0C6B", "#DC267F", "#FE6100", ...
    "#FFB000", "#016D24"];
% E_all(1,:) = E_all(1,:) - vs.alpha_i^2*vs.kD*vs.d_vals.^2/2;
% E_all(4,:) = E_all(4,:) - vs.alpha_i^2*vs.kD*vs.d_vals.^2/2;
for ii=1:6
%     plot(sigma_vals, E_all(ii,:)./(vs.kD.*vs.d_vals.^2), ...
%         strcat(colours(ii),lines(ii)))
%     plot(axes1,sigma_vals, E_all(ii,:), ...
%         strcat(colours(ii),lines(ii)))
    plot(axes1,sigma_vals, E_all(ii,:).*N_vals, ...
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
xlim([sigma_vals(1), sigma_vals(end)])
% ylim([min_E, max_E])

% plot second axes
axes2 = axes('Position',axes1.Position,'XAxisLocation','top','YAxisLocation','right','color','none')
hold(axes2, 'on');
% plot(axes2,N_vals,zeros(size(N_vals)),'-k')
set(axes2,'XTickLabel',round(N_vals(round(linspace(1,length(N_vals),5))),-5));
set(axes2,'ytick',linspace(0,1,5))
set(axes2,'ytick',[])
axes2.Color = 'none';
axes2.Box = 'off';
axes2.XLabel.String = '$N$';
axes2.XAxisLocation = 'top';
end