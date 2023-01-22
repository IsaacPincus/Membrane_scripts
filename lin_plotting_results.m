clear variables

tic
if isfile('combined.mat')
    load('combined.mat')
else
    num_tasks = 48;
    solution_set_total = {};
    E_all_total = [];
    alpha_A_vals_total = [];
    alpha_B_vals_total = [];
    phi_vals_total = [];
    h_phi_vals_total = [];
    S_A_vals_total = [];
    S_B_vals_total = [];
    Sigma_vals_total = [];
    for ii=0:num_tasks-1
        load_name = sprintf('data/task%i_results.mat', ii );
        data = load(load_name);
        solution_set_total = [solution_set_total, data.solution_set];
        E_all_total = cat(2,E_all_total, data.E_all);
        alpha_A_vals_total = horzcat(alpha_A_vals_total, data.alpha_A_vals);
        alpha_B_vals_total = horzcat(alpha_B_vals_total, data.alpha_B_vals);
        phi_vals_total = horzcat(phi_vals_total, data.phi_vals);
        h_phi_vals_total = horzcat(h_phi_vals_total, data.h_phi_vals);
        S_A_vals_total = horzcat(S_A_vals_total, data.S_A_vals);
        S_B_vals_total = horzcat(S_B_vals_total, data.S_B_vals);
        Sigma_vals_total = horzcat(Sigma_vals_total, data.Sigma_vals);
    end
    
    parameter_set = data.parameter_set;
    save('combined.mat')
end
toc

%% base case
% generate list of independent variables to run, which should be in the
% order [epsilon, n0, d, R, kD, kappa, alpha_i] for each row
R_vals = 1;                 					% um
sigma_vals = logspace(-3,-1.8,48);                % surface fraction
%d_vals = sqrt(R_vals.^2./sigma_vals);           % um
kD_vals = 300/10^12*1e9;                        % picoJ/um^2
kD_base = 300/10^12*1e9;                        % picoJ/um^2
zeta_vals = [0.01,0.02,0.03];                  % dimensionless
epsilon_vals = -zeta_vals*kD_base;              % picoJ/um^2
n0_vals = 1;                                    % fraction
kappa_vals = 1e-19*1e12;         % picoJ
% kappa_vals = logspace(-21,-15, 60)*1e12;
alpha_i_vals = [0.013,0.0197];                            % fraction
% other constants
ii = 0;
slice = [];
for rr = 1:length(R_vals)
    for ss = 1:length(sigma_vals)
        for kk = 1:length(kD_vals)
            for ee = 1:length(epsilon_vals)
                for nn = 1:length(n0_vals)
                    for pp = 1:length(kappa_vals)
                        for aa = 1:length(alpha_i_vals)
                            ii = ii+1;
                            if aa==1 && ee==2
                                slice = [slice, ii];
                            end
                        end
                    end
                end
            end
        end
    end
end
size_params = size(parameter_set);
base_index = slice(1);
% step = 1;
% if using R
% slice = base_index:step:size_params(1);
% if using alpha_i
% slice = (base_index-1)*step+1:base_index*step;
% base_index = (base_index-1)*step+1;
epsilon = parameter_set(base_index,1);
n0 = parameter_set(base_index,2);
d = parameter_set(base_index,3);
R = parameter_set(base_index,4);
kD = parameter_set(base_index,5);
kappa = parameter_set(base_index,6);
alpha_i = parameter_set(base_index,7);
zeta = epsilon*n0/kD;
sigma = R^2/d^2;

changing_value = 3;
% param_1_vals = R^2./parameter_set(slice,changing_value).^2;
param_1_vals = R.^2./parameter_set(slice,changing_value).^2;
E_all = E_all_total(:, slice);
alpha_A_vals = alpha_A_vals_total(slice);
alpha_B_vals = alpha_B_vals_total(slice);
phi_vals = phi_vals_total(slice);
h_phi_vals = h_phi_vals_total(slice);
S_A_vals = S_A_vals_total(slice);
S_B_vals = S_B_vals_total(slice);
Sigma_vals = Sigma_vals_total(slice);

solution_set = solution_set_total(slice);

% specific curves to plot
% plot_curves = [1,5, 10, 15, 25, 30, 35];
plot_curves = round(linspace(1,48,7));
% plot_curves = [1,54,74,76,82,85,88];

small_param = 0.005;
large_param = 0.015;

colours_cb = ["#27E0D8", "#648FFF", "#1D0C6B", "#DC267F", "#FE6100", ...
    "#FFB000", "#016D24"];

%% plotting
close all
init_vals = param_1_vals(1:40);
on = ones(size(init_vals));

% xaxis_name = '$\alpha_i$';
xaxis_name = '$\sigma$';
% xaxis_name = '$R$ ($\mu$m)';
% xaxis_name = '$\kappa$';

% xscale = 'linear';
xscale = 'log';

figure('Position',[400,100,800,600]);
hold on
axes1 = gca;
axes1.XScale = xscale;
% xlabel('$\kappa$ (pJ)')
% xlabel('$\zeta$')
% xlabel('surface fraction $\sigma$')
xlabel(xaxis_name)
ylabel('E (pJ)')
lines = ["-", ":", ":", "--", ":", "--"];
colours = ['k', 'b', 'r', 'r', 'g', 'g'];
for ii=1:6
%     plot(param_1_vals, E_all(ii,:).*param_1_vals', ...
%         strcat(colours(ii),lines(ii)))
    plot(param_1_vals, E_all(ii,:), ...
        strcat(colours(ii),lines(ii)))
end
min_E = min(E_all, [], 'all')*1.1;
max_E = max(E_all, [], 'all')*1.1;
% small_param = param_1_vals(10);
% large_param = param_1_vals(30);
patch([small_param, small_param, large_param, large_param],...
      [min_E, max_E, max_E, min_E], 'm', 'FaceAlpha', 0.05, ...
      'EdgeColor', 'None', 'HandleVisibility', 'off')
% plot(init_vals, on*data.E_total_init, 'k-.', 'linewidth', 0.5);
% plot(init_vals, on*data.E_adhesion_init, 'b-.', 'linewidth', 0.5);
% plot(init_vals, on*data.E_stretch_B_init, 'r-.', 'linewidth', 0.5);
% plot(init_vals, on*-0.001466887916365, 'k-.', 'linewidth', 0.5);
% plot(init_vals, on*-0.002422090201862, 'b-.', 'linewidth', 0.5);
% plot(init_vals, on*9.516558874055476e-04, 'r-.', 'linewidth', 0.5);
plot_vert_lines(param_1_vals(plot_curves), colours_cb)
legend({'$E_\mathrm{total}$','$E_\mathrm{adhesion}$',...
    '$E_\mathrm{stretch,A}$','$E_\mathrm{stretch,B}$',...
    '$E_\mathrm{bend,A}$','$E_\mathrm{bend,B}$'}, 'Box','off',...
    'location', 'best')
annotation('textbox', [0.3,0.6,0.4,0.2], 'String',...
    [sprintf('$k_D = %0.2g$ pJ/$\\mu$m$^2$ \n', kD),...
    sprintf('$\\kappa = %0.2g$ pJ \n', kappa),...
    sprintf('$\\alpha_i = %0.2g$ \n', alpha_i),...
    sprintf('$d = %0.2g$ $\\mu$m \n', d),...
    sprintf('$\\epsilon n_0 = %0.2g k_\\mathrm{D}$ \n', zeta)],...
    'interpreter', 'latex','FontSize',16,...
    'FitBoxToText','on','LineStyle','none', 'Color','b')
xlim([param_1_vals(1), param_1_vals(end)])
ylim([min_E, max_E])

figure('Position',[200,100,1100,600])
for ii=1:6
%     plot(param_1_vals, E_all(ii,:).*param_1_vals', ...
%         strcat(colours(ii),lines(ii)))
    subplot(2,3,ii)
    plot(param_1_vals, E_all(ii,:), ...
        strcat(colours(ii),lines(ii)))
    axes1 = gca;
    axes1.XScale = xscale;
%     xlabel(xaxis_name, 'FontSize',14)
%     ylabel('E (pJ)', 'FontSize',14)
%     axes1.YScale = 'log';
end

% single component linear and log scale
figure('Position',[400,100,800,600]);
hold on
axes1 = gca;
axes1.XScale = xscale;
% xlabel('$\kappa$ (pJ)')
% xlabel('$\zeta$')
% xlabel('surface fraction $\sigma$')
xlabel(xaxis_name)
ylabel('E (pJ)')
lines = ["-", ":", ":", "--", ":", "--"];
colours = ['k', 'b', 'r', 'r', 'g', 'g'];
for ii=2
    plot(param_1_vals, E_all(ii,:), ...
        strcat(colours(ii),lines(ii)))
end
legend({'$E_\mathrm{adhesion}$'}, 'Box','off',...
    'location', 'best')
xlim([param_1_vals(1), param_1_vals(end)])
ylim([min_E, max_E])

figure('Position',[400,100,800,600]);
hold on
axes1 = gca;
axes1.XScale = xscale;
axes1.YScale = 'log';
% xlabel('$\kappa$ (pJ)')
% xlabel('$\zeta$')
% xlabel('surface fraction $\sigma$')
xlabel(xaxis_name)
ylabel('E (pJ)')
lines = ["-", ":", ":", "--", ":", "--"];
colours = ['k', 'b', 'r', 'r', 'g', 'g'];
for ii=2
    plot(param_1_vals, E_all(ii,:), ...
        strcat(colours(ii),lines(ii)))
end
legend({'$E_\mathrm{adhesion}$'}, 'Box','off',...
    'location', 'best')
xlim([param_1_vals(1), param_1_vals(end)])
ylim([min_E, max_E])

%%

% positive percentages
figure('Position',[400,100,800,600]);
hold on
axes1 = gca;
axes1.XScale = xscale;
% xlabel('$\kappa$ (pJ)')
% xlabel('$\zeta$')
% xlabel('surface fraction $\sigma$')
xlabel(xaxis_name)
ylabel('Fraction of Energy')
lines = ["-", ":", ":", "--", ":", "--"];
colours = ['k', 'b', 'r', 'r', 'g', 'g'];
E_pos = E_all(1,:)-E_all(2,:);
E_pos_2 = sum(E_all(3:end,:));
E_perc = E_all./E_pos;
for ii=3:6
    plot(param_1_vals, E_perc(ii,:), ...
        strcat(colours(ii),lines(ii)))
end
plot(param_1_vals, E_perc(3,:)+E_perc(4,:), 'r-');
plot(param_1_vals, E_perc(5,:)+E_perc(6,:), 'g-');
axes_ylims = get(gca, 'YLim');
patch([small_param, small_param, large_param, large_param],...
      [axes_ylims(1), axes_ylims(2), axes_ylims(2), axes_ylims(1)],...
      'm', 'FaceAlpha', 0.05, ...
      'EdgeColor', 'None', 'HandleVisibility', 'off')
plot_vert_lines(param_1_vals(plot_curves), colours_cb)
legend({'$E_\mathrm{stretch,A}$','$E_\mathrm{stretch,B}$',...
    '$E_\mathrm{bend,A}$','$E_\mathrm{bend,B}$',...
    '$E_\mathrm{stretch,Total}$','$E_\mathrm{bend,Total}$'},...
    'Box','off','location', 'best')

%% areas
figure('Position',[400,100,800,600]);
hold on
axes1 = gca;
axes1.XScale = xscale;
% xlim([0,90]);
% xlabel('$\kappa$ (pJ)')
% xlabel('$\zeta$')
% xlabel('surface fraction $\sigma$')
xlabel(xaxis_name)
% yyaxis left
ylabel('Area')
% plot(init_vals, on*S_B_init, '-.', 'linewidth', 0.5, 'HandleVisibility','off');
% plot(init_vals, on*(S_B_init+S_A_init), '-.', 'linewidth', 0.5, 'HandleVisibility','off');
plot(param_1_vals, S_B_vals, '--', 'displayname', '$S_B$')
r_phi_vals = parameter_set(slice,4).*sin(phi_vals');
plot(param_1_vals, d^2-r_phi_vals.^2*pi, 'k:', 'displayname', '$d^2-\pi r_\phi^2$')
plot(param_1_vals, S_A_vals+S_B_vals, '-', 'displayname', '$S_\mathrm{total}$')
% yyaxis right
ylabel('Area')
% plot(init_vals, on*S_A_init, '-.', 'linewidth', 0.5, 'HandleVisibility','off');
plot(param_1_vals, S_A_vals, '--', 'displayname', '$S_A$')
axes_ylims = get(gca, 'YLim');
patch([small_param, small_param, large_param, large_param],...
      [axes_ylims(1), axes_ylims(2), axes_ylims(2), axes_ylims(1)],...
      'm', 'FaceAlpha', 0.05, ...
      'EdgeColor', 'None', 'HandleVisibility', 'off')
plot_vert_lines(param_1_vals(plot_curves), colours_cb)
legend

% total lipids
figure('Position',[400,100,800,600]);
hold on
axes1 = gca;
axes1.XScale = xscale;
% xlim([0,90]);
% xlabel('$\kappa$ (pJ)')
% xlabel('$\zeta$')
% xlabel('surface fraction $\sigma$')
xlabel(xaxis_name)
% yyaxis left
ylabel('Lipids')
% plot(init_vals, on*S_B_init, '-.', 'linewidth', 0.5, 'HandleVisibility','off');
% plot(init_vals, on*(S_B_init+S_A_init), '-.', 'linewidth', 0.5, 'HandleVisibility','off');
plot(param_1_vals, S_B_vals./(1+alpha_B_vals), '--', 'displayname', '$S_B$')
r_phi_vals = parameter_set(slice,4)'.*sin(phi_vals);
no_bend_B_lip = (d^2-r_phi_vals.^2*pi)./(1+alpha_B_vals);
plot(param_1_vals, no_bend_B_lip, 'k:', 'displayname', '$d^2-\pi r_\phi^2$')
plot(param_1_vals, S_A_vals./(1+alpha_A_vals)+S_B_vals./(1+alpha_B_vals),...
    '-', 'displayname', '$S_\mathrm{total}$')
% yyaxis right
ylabel('Lipids')
% plot(init_vals, on*S_A_init, '-.', 'linewidth', 0.5, 'HandleVisibility','off');
plot(param_1_vals, S_A_vals./(1+alpha_A_vals), '--', 'displayname', '$S_A$')
axes_ylims = get(gca, 'YLim');
patch([small_param, small_param, large_param, large_param],...
      [axes_ylims(1), axes_ylims(2), axes_ylims(2), axes_ylims(1)],...
      'm', 'FaceAlpha', 0.05, ...
      'EdgeColor', 'None', 'HandleVisibility', 'off')
plot_vert_lines(param_1_vals(plot_curves), colours_cb)
legend

%% individual total lipids
figure('Position',[400,100,800,600]);
hold on
axes1 = gca;
axes1.XScale = xscale;
axes1.YScale = 'log';
% xlim([0,90]);
% xlabel('$\kappa$ (pJ)')
% xlabel('$\zeta$')
% xlabel('surface fraction $\sigma$')
xlabel(xaxis_name)
% yyaxis left
ylabel('Lipids')
% plot(init_vals, on*S_B_init, '-.', 'linewidth', 0.5, 'HandleVisibility','off');
% plot(init_vals, on*(S_B_init+S_A_init), '-.', 'linewidth', 0.5, 'HandleVisibility','off');
% plot(param_1_vals, S_B_vals./(1+alpha_B_vals), '--', 'displayname', '$S_B$')
% r_phi_vals = parameter_set(slice,4)'.*sin(phi_vals);
% no_bend_B_lip = (d^2-r_phi_vals.^2*pi)./(1+alpha_B_vals);
% plot(param_1_vals, no_bend_B_lip, 'k:', 'displayname', '$d^2-\pi r_\phi^2$')
% plot(param_1_vals, S_A_vals./(1+alpha_A_vals)+S_B_vals./(1+alpha_B_vals),...
%     '-', 'displayname', '$S_\mathrm{total}$')
% yyaxis right
ylabel('Lipids')
% plot(init_vals, on*S_A_init, '-.', 'linewidth', 0.5, 'HandleVisibility','off');
plot(param_1_vals, S_A_vals./(1+alpha_A_vals), '--', 'displayname', '$S_A$')
axes_ylims = get(gca, 'YLim');
patch([small_param, small_param, large_param, large_param],...
      [axes_ylims(1), axes_ylims(2), axes_ylims(2), axes_ylims(1)],...
      'm', 'FaceAlpha', 0.05, ...
      'EdgeColor', 'None', 'HandleVisibility', 'off')
plot_vert_lines(param_1_vals(plot_curves), colours_cb)
legend

%% stretches
figure('Position',[400,100,800,600]);
hold on
axes1 = gca;
axes1.XScale = xscale;
% xlim([0,90]);
% xlabel('$\kappa$ (pJ)')
% xlabel('$\zeta$')
% xlabel('surface fraction $\sigma$')
xlabel(xaxis_name)
ylabel('$\alpha$')
plot(param_1_vals, alpha_B_vals, 'displayname', '$\alpha_B$')
plot(param_1_vals, alpha_A_vals, 'displayname', '$\alpha_A$')
plot(param_1_vals, zeros(size(param_1_vals)), 'k--', 'HandleVisibility', 'off')
% plot(init_vals, on*alpha_A_init, 'r-.', 'linewidth', 0.5, 'HandleVisibility','off');
% plot(init_vals, on*alpha_B_init, 'b-.', 'linewidth', 0.5, 'HandleVisibility','off');
axes_ylims = get(gca, 'YLim');
patch([small_param, small_param, large_param, large_param],...
      [axes_ylims(1), axes_ylims(2), axes_ylims(2), axes_ylims(1)],...
      'm', 'FaceAlpha', 0.05, ...
      'EdgeColor', 'None', 'HandleVisibility', 'off')
plot_vert_lines(param_1_vals(plot_curves), colours_cb)
legend

figure();
hold on
xlabel('concentration ($\mu$m/ml)')
ylabel('Tension (mN/m)')
concentration = 600*param_1_vals/0.016;
plot(concentration, alpha_B_vals*kD*1e3)

% phi and h_phi
figure('Position',[400,100,800,600]);
hold on
axes1 = gca;
axes1.XScale = xscale;
% xlim([0,90]);
% xlabel('$\kappa$ (pJ)')
% xlabel('$\zeta$')
% xlabel('surface fraction $\sigma$')
xlabel(xaxis_name)
yyaxis left
ylabel('$\phi$')
plot(param_1_vals, rad2deg(phi_vals), 'displayname', '$\phi$')
% plot(init_vals, on*rad2deg(phi_init), '-.', 'linewidth', 0.5, 'HandleVisibility','off')
yyaxis right
ylabel('$h_\phi$')
plot(param_1_vals, h_phi_vals, 'displayname', '$h_\phi$')
axes_ylims = get(gca, 'YLim');
patch([small_param, small_param, large_param, large_param],...
      [axes_ylims(1), axes_ylims(2), axes_ylims(2), axes_ylims(1)],...
      'm', 'FaceAlpha', 0.05, ...
      'EdgeColor', 'None', 'HandleVisibility', 'off')
plot_vert_lines(param_1_vals(plot_curves), colours_cb)
legend

% energy per lipid
figure('Position',[400,100,800,600]);
hold on
axes1 = gca;
axes1.XScale = xscale;
axes1.YScale = 'log';
% ylim([-1e-2, 1e-2])
% xlabel('$\kappa$ (pJ)')
% xlabel('$\zeta$')
% xlabel('surface fraction $\sigma$')
xlabel(xaxis_name)
ylabel('Energy per lipid')
lines = ["-", ":", ":", "--", ":", "--"];
colours = ['k', 'b', 'r', 'r', 'g', 'g'];
E_per_lipid(2,:) = -E_all(2,:)./S_A_vals.*(1+alpha_A_vals);
E_per_lipid([3,5],:) = E_all([3,5],:)./S_A_vals.*(1+alpha_A_vals);
E_per_lipid([4,6],:) = E_all([4,6],:)./S_B_vals.*(1+alpha_B_vals);
for ii=3:6
    plot(param_1_vals, E_per_lipid(ii,:), ...
        strcat(colours(ii),lines(ii)))
end
axes_ylims = get(gca, 'YLim');
patch([small_param, small_param, large_param, large_param],...
      [axes_ylims(1), axes_ylims(2), axes_ylims(2), axes_ylims(1)],...
      'm', 'FaceAlpha', 0.05, ...
      'EdgeColor', 'None', 'HandleVisibility', 'off')
plot_vert_lines(param_1_vals(plot_curves), colours_cb)
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

%% replotting each function

figure('Position',[400,100,800,600]);
hold on
axis equal
xlabel('$r$')
ylabel('$h$')
ylabel('$\int r \nabla^2 h dr$')
axes1 = gca;
% axes1.YScale = 'log';
N = 1e4;
counter = 0;
for ii = plot_curves
    counter = counter+1;
    index = slice(ii);
    epsilon = parameter_set(index,1);
    n0 = parameter_set(index,2);
    d = parameter_set(index,3);
    R = parameter_set(index,4);
    kD = parameter_set(index,5);
    kappa = parameter_set(index,6);
    alpha_i = parameter_set(index,7);
    zeta = epsilon*n0/kD;
    sigma = R^2/d^2;

    alpha_A = alpha_A_vals(ii);
    alpha_B = alpha_B_vals(ii);
    phi = phi_vals(ii);
    h_phi = h_phi_vals(ii);
    
    % get the shape of the free region
    Sigma = kD*alpha_B;
    lambda = sqrt(kappa/Sigma);
    r_phi = sin(phi)*R;
    r = linspace(r_phi, d/2,N);
    
    [h,C,S_B, ~, lap_h, hderiv] = ...
        free_shape_linear_fixed_h(r, r_phi, d, phi, kappa, Sigma, h_phi);

    C3_bar = C(3);
    C4_bar = C(4);

    first_term(:,counter) = exp((r-d/2)/lambda).*C3_bar/lambda^2.*besseli(0,r/lambda,1);
    second_term(:,counter) = exp((r_phi-r)/lambda).*C4_bar/lambda^2.*besselk(0,r/lambda,1);

    first_term_int(counter) = trapz(r,r.*first_term(:,counter)'.^2);
    second_term_int(counter) = trapz(r,r.*second_term(:,counter)'.^2);
    
    S_A = 2*pi*R^2*(1-cos(phi));
    E_adhesion = epsilon*n0*S_A./(1+alpha_A);
    E_stretch_A = kD/2*(alpha_A.^2*S_A./(1+alpha_A));
    E_stretch_B = kD/2*(alpha_B.^2*S_B./(1+alpha_B));
    E_bend_B = kappa/2*2*pi*trapz(r, r.*lap_h.^2);
    E_bend_A = 4*pi*kappa*(1-cos(phi));
    E = E_adhesion + E_stretch_A + E_stretch_B + E_bend_A + E_bend_B;
    
    z = 1-cos(phi);
    k = sin(phi);
    h_lim = R*(z-1+lambda/R*k/(1-z)*...
        (besselk(0,r_phi/lambda)-besselk(0,r/lambda))/besselk(1,r_phi/lambda));

    h_lim = h_lim - (h_lim(1)-h(1));

%     h1 = plot(r, h, 'displayname', sprintf('$k_D/k_{D,0} = %0.2e$', ...
%         kD/(300/10^12*1e9)));
%     h1 = plot(r, h, 'displayname', sprintf('$\\alpha_i = %0.2e$', ...
%         alpha_i), 'Color', colours_cb(counter));
    h1 = plot(r, h, 'displayname', sprintf('$\\kappa = %0.2e$', ...
        kappa), 'Color', colours_cb(counter));
%     h2 = plot(r, h_lim, '--', 'displayname', sprintf('$R = %0.2e$', ...
%         R), 'Color', colours_cb(counter), 'HandleVisibility','off');
%     yyaxis left
%     h1 = plot(r, r.*lap_h.^2,'-', 'displayname', sprintf('$R = %0.2e$', ...
%         R), 'Color', colours_cb(counter));
%     yyaxis right
%     h2 = plot(r, cumtrapz(r,r.*lap_h.^2),'-',  'displayname', sprintf('$R = %0.2e$', ...
%         R), 'Color', colours_cb(counter));
    colour = h1.Color;
    t = linspace(-pi/2,-pi/2+phi,1000);
    x = cos(t)*R;
    % y = sin(t)*R+(R*cos(phi)+h(1));
    y = sin(t)*R+R*cos(phi)+h(1);
    plot(x,y, ':', 'HandleVisibility','off', 'color', colours_cb(counter))
end

legend

%% plotting solution space
% figure();
% hold on
% ylim([min(cat(1,solutions.Fval))*1.1, 1e-5])
% plot(arrayfun(@(x)x.Fval,solutions),'k*')
% xlabel('Solution number')
% ylabel('Function value')
% % title('Solution Function Values')

solutions = solution_set{1};
out = solutions(1).X;
fval2 = solutions(1).Fval;

% zlim([-1e-3, 1e-5])
for ii = 1:size(solutions,2)
    X0 = solutions(ii).X0{1};
    X = solutions(ii).X;
    alpha_A_val(ii) = X(1);
    alpha_B_val(ii) = X(2);
    phi_val(ii) = X(3);
    h_phi_val(ii) = X(4);
    E_val(ii) = solutions(ii).Fval;
    dist(ii) = vecnorm(out-X);
    E_dist(ii) = E_val(ii) - fval2;
%     plot3(rad2deg(phi_val(ii)), h_phi_val(ii), E_val(ii), 'k.')
end
min_energy = min(cat(1,solutions.Fval))*1.1;
max_energy = 2e-4;
fig = figure('Position',[300,75,900,700]);
s_plot = subplot(2,2,1);
hold on
xlim([0,90])
ylim([min_energy, max_energy])
% plot(rad2deg(phi_val), E_val, 'o');
scatter(rad2deg(phi_val), E_val,...
    arrayfun(@(x)log(length(x.X0))*20+50,solutions),...
    arrayfun(@(x)length(x.X0),solutions), 'filled',...
    'Marker', 'o');
% c_bar = colorbar;
% c_bar.TickLabelInterpreter = 'latex';
% c_bar.Label.String = 'Multiplicity';
% colormap('winter')
% set(gca,'ColorScale','log')
% surf(phi_val, h_phi_val, E_val);
% plot(arrayfun(@(x)x.Fval,solutions),'k*')
xlabel('$\phi$')
ylabel('$E$')
pos = get(s_plot, 'Position'); 
posnew = pos;
posnew(1) = posnew(1) - 0.03;
set(s_plot, 'Position', posnew);

% figure();
s_plot = subplot(2,2,2);
hold on
xlim([-2*R, 2*R])
ylim([min_energy, max_energy])
% plot(rad2deg(phi_val), E_val, 'o');
scatter(h_phi_val, E_val,...
    arrayfun(@(x)log(length(x.X0))*20+50,solutions),...
    arrayfun(@(x)length(x.X0),solutions), 'filled',...
    'Marker', 'o');
% c_bar = colorbar;
% c_bar.TickLabelInterpreter = 'latex';
% c_bar.Label.String = 'Multiplicity';
% colormap('winter')
% set(gca,'ColorScale','log')
% surf(phi_val, h_phi_val, E_val);
% plot(arrayfun(@(x)x.Fval,solutions),'k*')
xlabel('$h_\phi$')
ylabel('$E$')
pos = get(s_plot, 'Position'); 
posnew = pos;
posnew(1) = posnew(1) - 0.03;
set(s_plot, 'Position', posnew);

% figure();
s_plot = subplot(2,2,3);
hold on
xlim([-0.1, 0.1])
ylim([min_energy, max_energy])
% plot(rad2deg(phi_val), E_val, 'o');
scatter(alpha_A_val, E_val,...
    arrayfun(@(x)log(length(x.X0))*20+50,solutions),...
    arrayfun(@(x)length(x.X0),solutions), 'filled',...
    'Marker', 'o');
% c_bar = colorbar;
% c_bar.TickLabelInterpreter = 'latex';
% c_bar.Label.String = 'Multiplicity';
% colormap('winter')
% set(gca,'ColorScale','log')
% surf(phi_val, h_phi_val, E_val);
% plot(arrayfun(@(x)x.Fval,solutions),'k*')
xlabel('$\alpha_\mathrm{A}$')
ylabel('$E$')
pos = get(s_plot, 'Position'); 
posnew = pos;
posnew(1) = posnew(1) - 0.03;
set(s_plot, 'Position', posnew);

% figure();
s_plot = subplot(2,2,4);
hold on
xlim([-0.1, 0.1])
ylim([min_energy, max_energy])
% plot(rad2deg(phi_val), E_val, 'o');
scatter(alpha_B_val, E_val,...
    arrayfun(@(x)log(length(x.X0))*20+50,solutions),...
    arrayfun(@(x)length(x.X0),solutions), 'filled',...
    'Marker', 'o');
% c_bar = colorbar;
% c_bar.TickLabelInterpreter = 'latex';
% c_bar.Label.String = 'Multiplicity';
% colormap('winter')
% set(gca,'ColorScale','log')
% surf(phi_val, h_phi_val, E_val);
% plot(arrayfun(@(x)x.Fval,solutions),'k*')
xlabel('$\alpha_\mathrm{B}$')
ylabel('$E$')
pos = get(s_plot, 'Position'); 
posnew = pos;
posnew(1) = posnew(1) - 0.03;
set(s_plot, 'Position', posnew);

h = axes(fig,'visible','off'); 
h.Toolbar.Visible = 'off';
c_bar = colorbar(h,'Position',[0.90 0.168 0.022 0.7]);  % attach colorbar to h
set(h,'ColorScale','log')
c_bar.TickLabelInterpreter = 'latex';
c_bar.Label.String = 'Multiplicity';
colormap('winter')
c_bar.Label.Position = [1.5,0.35];
set(fig, 'Children', flipud(fig.Children))
% colormap(c,'jet')
% caxis(h,[minColorLimit,maxColorLimit]);             % set colorbar limits

% figure();
% hold on
% xlim([0,90])
% ylim([-R, R])
% zlim([min_energy, max_energy])
% scatter3(rad2deg(phi_val), h_phi_val, E_val,...
%     arrayfun(@(x)log(length(x.X0))*20+50,solutions),...
%     arrayfun(@(x)length(x.X0),solutions), 'filled',...
%     'Marker', 'o');
% c_bar = colorbar;
% c_bar.TickLabelInterpreter = 'latex';
% c_bar.Label.String = 'Multiplicity';
% colormap('winter')
% set(gca,'ColorScale','log')
% xlabel('$\phi$')
% ylabel('$h_\phi$')
% zlabel('$E$')

function plot_vert_lines(x_vals, colours)
    for ii=1:length(x_vals)
        plot(ones(1,2)*x_vals(ii), [get(gca, 'YLim')],...
            ':o', 'linewidth', 0.5, 'HandleVisibility','off',...
            'Color',colours(ii), 'MarkerFaceColor',colours(ii));
    end
end