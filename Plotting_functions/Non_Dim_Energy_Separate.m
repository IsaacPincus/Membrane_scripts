function fig_handle = Non_Dim_Energy_Separate(xaxis_name, xscale, yaxis_name, yscale,...
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
symbols = ['p', 'o', 'o', 's', 'o', 's'];
colours_cb = ["#27E0D8", "#648FFF", "#1D0C6B", "#DC267F", "#FE6100", ...
    "#FFB000", "#016D24"];
markersize = 2;
E_all_ana = [vs.E_all(1,:)./(vs.d_vals.^2*vs.kD);...
    vs.zeta./(1+vs.alpha_A_vals)*2*pi.*(1-cos(vs.phi_vals)).*vs.sigma_vals;...
    0.5*vs.alpha_A_vals.^2./(1+vs.alpha_A_vals)*2*pi.*(1-cos(vs.phi_vals)).*vs.sigma_vals;...
    0.5*vs.alpha_B_vals.^2./(1+vs.alpha_B_vals).*(1-pi*vs.sigma_vals.*sin(vs.phi_vals).^2);...
    4*pi*vs.kappa_vals./(vs.d_vals.^2*vs.kD).*(1-cos(vs.phi_vals));...
    pi/2*vs.alpha_B_vals.*vs.sigma_vals.*tan(vs.phi_vals).^2.*sin(vs.phi_vals).^2.*...
    (besselk(1,vs.R_vals./vs.lambda_vals.*sin(vs.phi_vals),1).^2-...
     besselk(0,vs.R_vals./vs.lambda_vals.*sin(vs.phi_vals),1).^2)./...
     besselk(1,vs.R_vals./vs.lambda_vals.*sin(vs.phi_vals),1).^2];

names = {'$E_\mathrm{total}$','$E_\mathrm{adhesion}$',...
    '$E_\mathrm{stretch,A}$','$E_\mathrm{stretch,B}$',...
    '$E_\mathrm{bend,A}$','$E_\mathrm{bend,B}$'};
for ii=1:6
    subplot(2,3,ii)
    hold on
    plot(vs.param_1_vals, E_all_ana(ii,:), ...
        strcat('y',symbols(ii)),'markersize', markersize)
    plot(vs.param_1_vals, vs.E_all(ii,:)./(vs.d_vals.^2*vs.kD), ...
        strcat(colours(ii),lines(ii)))
    axes1 = gca;
    axes1.XScale = xscale;
    xlabel(xaxis_name)
    ylabel(yaxis_name)
    title(names{ii})
    axes1.YScale = yscale;
end
xlim([vs.param_1_vals(1), vs.param_1_vals(end)])
end