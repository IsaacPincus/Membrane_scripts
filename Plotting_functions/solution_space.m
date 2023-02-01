function fig_handle = solution_space(xaxis_name, xscale, yaxis_name, yscale,...
    vs, patch_lims, solution_set, solution_number)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% fig_handle = figure('Position',[400,100,800,600]);
% figure();
% hold on
% ylim([min(cat(1,solutions.Fval))*1.1, 1e-5])
% plot(arrayfun(@(x)x.Fval,solutions),'k*')
% xlabel('Solution number')
% ylabel('Function value')
% % title('Solution Function Values')

solutions = solution_set{solution_number};
out = solutions(1).X;
fval2 = solutions(1).Fval;
R = vs.R_vals(solution_number);

alpha_A_val = zeros(size(solutions));
alpha_B_val = zeros(size(solutions));
phi_val = zeros(size(solutions));
h_phi_val = zeros(size(solutions));
E_val = zeros(size(solutions));
dist = zeros(size(solutions));
E_dist = zeros(size(solutions));

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
min_energy = min(cat(1,solutions.Fval));
if min_energy<0
    min_energy = min_energy*1.1;
else
    min_energy = min_energy*0.9;
end
max_energy = 2e-4;
fig_handle = figure('Position',[300,75,900,700]);
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

h = axes(fig_handle,'visible','off'); 
h.Toolbar.Visible = 'off';
c_bar = colorbar(h,'Position',[0.90 0.168 0.022 0.7]);  % attach colorbar to h
set(h,'ColorScale','log')
c_bar.TickLabelInterpreter = 'latex';
c_bar.Label.String = 'Multiplicity';
colormap('winter')
c_bar.Label.Position = [1.5,0.35];
set(fig_handle, 'Children', flipud(fig_handle.Children))
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

end