function fig_handle = stretches_ratio(xaxis_name, xscale, yaxis_name, yscale,...
    vs, patch_lims, plot_curves, plot_Baulin, fig_handle_inp)
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

% plot(vs.param_1_vals, vs.alpha_B_vals./vs.alpha_A_vals, 'displayname', '$\alpha_B/\alpha_A$')
plot(vs.param_1_vals, (1+vs.alpha_B_vals)./(1+vs.alpha_A_vals), 'displayname', '$(1+\alpha_B)/(1+\alpha_A$)')
% plot(vs.param_1_vals, , 'displayname', '$\alpha_A$')



legend

if ~exist("plot_Baulin", 'var') || plot_Baulin==0
    plot_Baulin = 0;
else
    opts = detectImportOptions("Baulin_Fig1C.csv");
    data = readmatrix('Baulin_Fig1C.csv', opts);
%     opts = detectImportOptions("Baulin_Fig1D.csv");
%     data = readmatrix('Baulin_Fig1D.csv', opts);
    % X = data(:,1:2:end);
    % Y = data(:,2:2:end);
    % plot(X, Y, 'displayName', 'Baulin 2021')
    X = data(:,1:2:end);
    Y = data(:,2:2:end);
    plot(X, Y, 'displayName', 'Baulin 2021, $\zeta = -0.005$')
end

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