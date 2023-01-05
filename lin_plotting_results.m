clear variables

num_tasks = 48;
solution_set = {};
E_all = [];
alpha_A_vals = [];
alpha_B_vals = [];
phi_vals = [];
h_phi_vals = [];
S_A_vals = [];
S_B_vals = [];
Sigma_vals = [];
for ii=0:num_tasks-1
    load_name = sprintf('data/task%i_results.mat', ii);
    data = load(load_name);
    solution_set = [solution_set, data.solution_set];
    E_all = cat(2,E_all, data.E_all);
    alpha_A_vals = horzcat(alpha_A_vals, data.alpha_A_vals);
    alpha_B_vals = horzcat(alpha_B_vals, data.alpha_B_vals);
    phi_vals = horzcat(phi_vals, data.phi_vals);
    h_phi_vals = horzcat(h_phi_vals, data.h_phi_vals);
    S_A_vals = horzcat(S_A_vals, data.S_A_vals);
    S_B_vals = horzcat(S_B_vals, data.S_B_vals);
    Sigma_vals = horzcat(Sigma_vals, data.Sigma_vals);
end

parameter_set = data.parameter_set;

% base case
base_index = 1;
epsilon = parameter_set(base_index,1);
n0 = parameter_set(base_index,2);
d = parameter_set(base_index,3);
R = parameter_set(base_index,4);
kD = parameter_set(base_index,5);
kappa = parameter_set(base_index,6);
alpha_i = parameter_set(base_index,7);
zeta = epsilon*n0/kD;
sigma = sqrt(R^2/d^2);
%
changing_value = 6;
param_1_vals = parameter_set(:,changing_value)/1e-7;

%% plotting
close all
init_vals = param_1_vals(1:40);
on = ones(size(init_vals));

xaxis_name = '$\kappa/\kappa_0$';

figure('Position',[400,100,800,600]);
hold on
axes1 = gca;
axes1.XScale = 'log';
% xlabel('$\kappa$ (pJ)')
% xlabel('$\zeta$')
% xlabel('surface fraction $\sigma$')
xlabel(xaxis_name)
ylabel('E (pJ)')
lines = ["-", ":", ":", "--", ":", "--"];
colours = ['k', 'b', 'r', 'r', 'g', 'g'];
for ii=1:6
    plot(param_1_vals, E_all(ii,:), ...
        strcat(colours(ii),lines(ii)))
end
% plot(init_vals, on*data.E_total_init, 'k-.', 'linewidth', 0.5);
% plot(init_vals, on*data.E_adhesion_init, 'b-.', 'linewidth', 0.5);
% plot(init_vals, on*data.E_stretch_B_init, 'r-.', 'linewidth', 0.5);
plot(init_vals, on*-0.001466887916365, 'k-.', 'linewidth', 0.5);
plot(init_vals, on*-0.002422090201862, 'b-.', 'linewidth', 0.5);
plot(init_vals, on*9.516558874055476e-04, 'r-.', 'linewidth', 0.5);
legend({'$E_\mathrm{total}$','$E_\mathrm{adhesion}$',...
    '$E_\mathrm{stretch,A}$','$E_\mathrm{stretch,B}$',...
    '$E_\mathrm{bend,A}$','$E_\mathrm{bend,B}$'}, 'Box','off',...
    'location', 'best')
annotation('textbox', [0.3,0.6,0.4,0.2], 'String',...
    [sprintf('$R = %0.2g$ $\\mu$m \n', R),...
    sprintf('$k_D = %0.2g$ pJ/$\\mu$m$^2$ \n', kD),...
    sprintf('$\\alpha_i = %0.2g$ \n', alpha_i),...
    sprintf('$d^2 = %0.2g R^2$ \n', sigma),...
    sprintf('$\\epsilon n_0 = %0.2g k_\\mathrm{D}$ \n', zeta)],...
    'interpreter', 'latex','FontSize',16,...
    'FitBoxToText','on','LineStyle','none', 'Color','b')

% positive percentages
figure('Position',[400,100,800,600]);
hold on
axes1 = gca;
axes1.XScale = 'log';
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
legend({'$E_\mathrm{stretch,A}$','$E_\mathrm{stretch,B}$',...
    '$E_\mathrm{bend,A}$','$E_\mathrm{bend,B}$',...
    '$E_\mathrm{stretch,Total}$','$E_\mathrm{bend,Total}$'},...
    'Box','off','location', 'best')

% areas
figure('Position',[400,100,800,600]);
hold on
axes1 = gca;
axes1.XScale = 'log';
% xlim([0,90]);
% xlabel('$\kappa$ (pJ)')
% xlabel('$\zeta$')
% xlabel('surface fraction $\sigma$')
xlabel(xaxis_name)
yyaxis left
ylabel('Area')
% plot(init_vals, on*S_B_init, '-.', 'linewidth', 0.5, 'HandleVisibility','off');
% plot(init_vals, on*(S_B_init+S_A_init), '-.', 'linewidth', 0.5, 'HandleVisibility','off');
plot(param_1_vals, S_B_vals, '--', 'displayname', '$S_B$')
plot(param_1_vals, S_A_vals+S_B_vals, '-', 'displayname', '$S_\mathrm{total}$')
yyaxis right
ylabel('Area')
% plot(init_vals, on*S_A_init, '-.', 'linewidth', 0.5, 'HandleVisibility','off');
plot(param_1_vals, S_A_vals, '--', 'displayname', '$S_A$')
legend

% stretches
figure('Position',[400,100,800,600]);
hold on
axes1 = gca;
axes1.XScale = 'log';
% xlim([0,90]);
% xlabel('$\kappa$ (pJ)')
% xlabel('$\zeta$')
% xlabel('surface fraction $\sigma$')
xlabel(xaxis_name)
ylabel('$\alpha$')
plot(param_1_vals, alpha_B_vals, 'displayname', '$\alpha_B$')
plot(param_1_vals, alpha_A_vals, 'displayname', '$\alpha_A$')
% plot(init_vals, on*alpha_A_init, 'r-.', 'linewidth', 0.5, 'HandleVisibility','off');
% plot(init_vals, on*alpha_B_init, 'b-.', 'linewidth', 0.5, 'HandleVisibility','off');
legend

% phi and h_phi
figure('Position',[400,100,800,600]);
hold on
axes1 = gca;
axes1.XScale = 'log';
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
legend

% energy per lipid
figure('Position',[400,100,800,600]);
hold on
axes1 = gca;
axes1.XScale = 'log';
axes1.YScale = 'log';
% ylim([-1e-2, 1e-2])
% xlabel('$\kappa$ (pJ)')
% xlabel('$\zeta$')
% xlabel('surface fraction $\sigma$')
xlabel(xaxis_name)
ylabel('Energy per lipid')
lines = ["-", ":", ":", "--", ":", "--"];
colours = ['k', 'b', 'r', 'r', 'g', 'g'];
E_per_lipid([2,3,5],:) = E_all([2,3,5],:)./S_A_vals.*(1+alpha_A_vals);
E_per_lipid([4,6],:) = E_all([4,6],:)./S_B_vals.*(1+alpha_B_vals);
for ii=3:6
    plot(param_1_vals, E_per_lipid(ii,:), ...
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

%% replotting each function

figure('Position',[400,100,800,600]);
hold on
axis equal
xlabel('$r$')
ylabel('$h$')
N = 1e6;
for ii = [1,20,40,60,80,96]
    epsilon = parameter_set(ii,1);
    n0 = parameter_set(ii,2);
    d = parameter_set(ii,3);
    R = parameter_set(ii,4);
    kD = parameter_set(ii,5);
    kappa = parameter_set(ii,6);
    alpha_i = parameter_set(ii,7);

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
    
    S_A = 2*pi*R^2*(1-cos(phi));
    E_adhesion = epsilon*n0*S_A./(1+alpha_A);
    E_stretch_A = kD/2*(alpha_A.^2*S_A./(1+alpha_A));
    E_stretch_B = kD/2*(alpha_B.^2*S_B./(1+alpha_B));
    E_bend_B = kappa/2*2*pi*trapz(r, r.*lap_h.^2);
    E_bend_A = 4*pi*kappa*(1-cos(phi));
    E = E_adhesion + E_stretch_A + E_stretch_B + E_bend_A + E_bend_B;

    h1 = plot(r, h, 'displayname', sprintf('$\\kappa/\\kappa_0 = %0.2e$', kappa/1e-7));
    colour = h1.Color;
    t = linspace(-pi/2,-pi/2+phi,1000);
    x = cos(t)*R;
    % y = sin(t)*R+(R*cos(phi)+h(1));
    y = sin(t)*R+R*cos(phi)+h(1);
    plot(x,y, ':', 'HandleVisibility','off', 'color', colour)
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
xlim([-R, R])
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
c_bar = colorbar(h,'Position',[0.90 0.168 0.022 0.7]);  % attach colorbar to h
set(h,'ColorScale','log')
c_bar.TickLabelInterpreter = 'latex';
c_bar.Label.String = 'Multiplicity';
colormap('winter')
c_bar.Label.Position = [1.5,0.35];
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