clear variables

load_name = 'testing.mat';
load(load_name);

%% plotting
init_vals = param_1_vals(1:20);
on = ones(size(init_vals));

figure('Position',[400,100,800,600]);
hold on
axes1 = gca;
axes1.XScale = 'log';
% xlabel('$\kappa$ (pJ)')
% xlabel('$\zeta$')
xlabel('surface fraction $\sigma$')
ylabel('E (pJ)')
lines = ["-", ":", ":", "--", ":", "--"];
colours = ['k', 'b', 'r', 'r', 'g', 'g'];
for ii=1:6
    plot(param_1_vals, E_all(ii,:), ...
        strcat(colours(ii),lines(ii)))
end
% plot(init_vals, on*E_total_init, 'k-.', 'linewidth', 0.5);
% plot(init_vals, on*E_adhesion_init, 'b-.', 'linewidth', 0.5);
% plot(init_vals, on*E_stretch_B_init, 'r-.', 'linewidth', 0.5);
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
xlabel('surface fraction $\sigma$')
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
legend({'$E_\mathrm{stretch,A}$','$E_\mathrm{stretch,B}$',...
    '$E_\mathrm{bend,A}$','$E_\mathrm{bend,B}$'}, 'Box','off',...
    'location', 'best')

% areas
figure('Position',[400,100,800,600]);
hold on
axes1 = gca;
axes1.XScale = 'log';
% xlim([0,90]);
% xlabel('$\kappa$ (pJ)')
% xlabel('$\zeta$')
xlabel('surface fraction $\sigma$')
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
xlabel('surface fraction $\sigma$')
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
xlabel('surface fraction $\sigma$')
yyaxis left
ylabel('$\phi$')
plot(param_1_vals, rad2deg(phi_vals), 'displayname', '$\phi$')
% plot(init_vals, on*rad2deg(phi_init), '-.', 'linewidth', 0.5, 'HandleVisibility','off')
yyaxis right
ylabel('$h_\phi$')
plot(param_1_vals, h_phi_vals, 'displayname', '$h_\phi$')
legend

%% replotting each function

figure('Position',[400,100,700,500]);
hold on
axis equal
xlabel('$r$')
ylabel('$h$')
for ii = [1,5,10,20,30,40]
    sigma = param_1_vals(ii);
    d = sqrt(R^2/sigma);    % um

    alpha_A = alpha_A_vals(ii);
    alpha_B = alpha_B_vals(ii);
    phi = phi_vals(ii);
    h_phi = h_phi_vals(ii);
    
    % get the shape of the free region
    Sigma = kD*alpha_B;
    lambda = sqrt(kappa/Sigma);
    r_phi = sin(phi)*R;
    r = linspace(r_phi, d/2,N);
    
    [h,C,S_B, ~, lap_h, hderiv] = free_shape_linear_fixed_h(r, r_phi, d, phi, kappa, Sigma, h_phi);

    h1 = plot(r, h, 'displayname', sprintf('$\\kappa = %0.2e$', kappa));
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
% figure();
subplot(2,2,1)
hold on
xlim([0,90])
ylim([min_energy, max_energy])
% plot(rad2deg(phi_val), E_val, 'o');
scatter(rad2deg(phi_val), E_val,...
    arrayfun(@(x)log(length(x.X0))*20+50,solutions),...
    arrayfun(@(x)length(x.X0),solutions), 'filled',...
    'Marker', 'o');
c_bar = colorbar;
c_bar.TickLabelInterpreter = 'latex';
c_bar.Label.String = 'Multiplicity';
colormap('winter')
set(gca,'ColorScale','log')
% surf(phi_val, h_phi_val, E_val);
% plot(arrayfun(@(x)x.Fval,solutions),'k*')
xlabel('$\phi$')
ylabel('$E$')

% figure();
subplot(2,2,2)
hold on
xlim([-R, R])
ylim([min_energy, max_energy])
% plot(rad2deg(phi_val), E_val, 'o');
scatter(h_phi_val, E_val,...
    arrayfun(@(x)log(length(x.X0))*20+50,solutions),...
    arrayfun(@(x)length(x.X0),solutions), 'filled',...
    'Marker', 'o');
c_bar = colorbar;
c_bar.TickLabelInterpreter = 'latex';
c_bar.Label.String = 'Multiplicity';
colormap('winter')
set(gca,'ColorScale','log')
% surf(phi_val, h_phi_val, E_val);
% plot(arrayfun(@(x)x.Fval,solutions),'k*')
xlabel('$h_\phi$')
ylabel('$E$')

% figure();
subplot(2,2,3)
hold on
xlim([-0.1, 0.1])
ylim([min_energy, max_energy])
% plot(rad2deg(phi_val), E_val, 'o');
scatter(alpha_A_val, E_val,...
    arrayfun(@(x)log(length(x.X0))*20+50,solutions),...
    arrayfun(@(x)length(x.X0),solutions), 'filled',...
    'Marker', 'o');
c_bar = colorbar;
c_bar.TickLabelInterpreter = 'latex';
c_bar.Label.String = 'Multiplicity';
colormap('winter')
set(gca,'ColorScale','log')
% surf(phi_val, h_phi_val, E_val);
% plot(arrayfun(@(x)x.Fval,solutions),'k*')
xlabel('$\alpha_\mathrm{A}$')
ylabel('$E$')

% figure();
subplot(2,2,4)
hold on
xlim([-0.1, 0.1])
ylim([min_energy, max_energy])
% plot(rad2deg(phi_val), E_val, 'o');
scatter(alpha_B_val, E_val,...
    arrayfun(@(x)log(length(x.X0))*20+50,solutions),...
    arrayfun(@(x)length(x.X0),solutions), 'filled',...
    'Marker', 'o');
c_bar = colorbar;
c_bar.TickLabelInterpreter = 'latex';
c_bar.Label.String = 'Multiplicity';
colormap('winter')
set(gca,'ColorScale','log')
% surf(phi_val, h_phi_val, E_val);
% plot(arrayfun(@(x)x.Fval,solutions),'k*')
xlabel('$\alpha_\mathrm{B}$')
ylabel('$E$')

figure();
hold on
xlim([0,90])
ylim([-R, R])
zlim([min_energy, max_energy])
scatter3(rad2deg(phi_val), h_phi_val, E_val,...
    arrayfun(@(x)log(length(x.X0))*20+50,solutions),...
    arrayfun(@(x)length(x.X0),solutions), 'filled',...
    'Marker', 'o');
c_bar = colorbar;
c_bar.TickLabelInterpreter = 'latex';
c_bar.Label.String = 'Multiplicity';
colormap('winter')
set(gca,'ColorScale','log')
xlabel('$\phi$')
ylabel('$h_\phi$')
zlabel('$E$')