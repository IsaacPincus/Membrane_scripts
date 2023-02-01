function fig_handle = Energy_Subplots(xaxis_name, xscale, yaxis_name, yscale,...
    vs, patch_lims, plot_curves, anno_string)
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
names = {'$E_\mathrm{total}$','$E_\mathrm{adhesion}$',...
    '$E_\mathrm{stretch,A}$','$E_\mathrm{stretch,B}$',...
    '$E_\mathrm{bend,A}$','$E_\mathrm{bend,B}$'};
for ii=1:6
    subplot(2,3,ii)
    plot(vs.param_1_vals, vs.E_all(ii,:), ...
        strcat(colours(ii),lines(ii)))
    axes1 = gca;
    axes1.XScale = xscale;
    axes1.YScale = yscale;
    axes1.XAxis.FontSize = 14;
    axes1.YAxis.FontSize = 14;
    xlabel(xaxis_name, 'FontSize',14)
    ylabel(yaxis_name, 'FontSize',14)
    title(names{ii}, 'FontSize',14)
end
xlim([vs.param_1_vals(1), vs.param_1_vals(end)])
end