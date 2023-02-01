function fig_handle = Non_Dim_Separated_Energies(xaxis_name, xscale, yaxis_name, yscale,...
    vs, patch_lims, plot_curves, curve_to_plot)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

fig_handle = figure('Position',[400,100,800,600]);

hold on
axes1 = gca;
axes1.XScale = xscale;
axes1.YScale = yscale;
xlabel(xaxis_name)
ylabel(yaxis_name)
range = 1:length(vs.param_1_vals);
colours_cb = ["#27E0D8", "#648FFF", "#1D0C6B", "#DC267F", "#FE6100", ...
    "#FFB000", "#016D24"];

if curve_to_plot==2
    % for adhesion energy
    plot(vs.param_1_vals(range), -vs.E_all(2,range)./(vs.d_vals(range).^2*vs.kD),...
        'k-','displayname', '$-E_\mathrm{adhesion}$');
    plot(vs.param_1_vals(range), 2*pi*(1-cos(vs.phi_vals(range))),...
        '-','displayname', '$2 \pi (1-\cos{\phi})$', 'color', colours_cb(1));
    plot(vs.param_1_vals(range), vs.sigma_vals(range),...
        'r-','displayname', '$\sigma$', 'color', colours_cb(4));
    plot(vs.param_1_vals(range), 1./(1+vs.alpha_A_vals(range)),...
        'g-','displayname', '$1/(1+\alpha_A)$', 'color', colours_cb(5));
    plot(vs.param_1_vals(range), -vs.zeta_vals(range),...
        'c-','displayname', '$-\zeta$', 'color', colours_cb(7));
    legend


elseif curve_to_plot==3
    % for stretching energy in A
    plot(vs.param_1_vals(range), vs.E_all(3,range)./(vs.d_vals(range).^2*vs.kD),...
        'k-','displayname', '$E_\mathrm{Stretch,A}$');
    plot(vs.param_1_vals(range), pi*(1-cos(vs.phi_vals(range))),...
        '-','displayname', '$\pi (1-\cos{\phi})$', 'color', colours_cb(1));
    plot(vs.param_1_vals(range), vs.sigma_vals(range),...
        'r-','displayname', '$\sigma$', 'color', colours_cb(4));
    plot(vs.param_1_vals(range), 1./(1+vs.alpha_A_vals(range)),...
        'g-','displayname', '$1/(1+\alpha_A)$', 'color', colours_cb(5));
    plot(vs.param_1_vals(range), vs.alpha_A_vals(range).^2,...
        'c-','displayname', '$\alpha_A^2$', 'color', colours_cb(7));
    legend

elseif curve_to_plot==4
    % for stretching energy in B
    plot(vs.param_1_vals(range), vs.E_all(4,range)./(vs.d_vals(range).^2*vs.kD),...
        'k-','displayname', '$E_\mathrm{Stretch,B}$');
    plot(vs.param_1_vals(range), pi*(1-cos(vs.phi_vals(range))),...
        '-','displayname', '$\pi (1-\cos{\phi})$', 'color', colours_cb(1));
    plot(vs.param_1_vals(range), vs.sigma_vals(range),...
        'r-','displayname', '$\sigma$', 'color', colours_cb(4));
    plot(vs.param_1_vals(range), 1./(1+vs.alpha_B_vals(range)),...
        'g-','displayname', '$1/(1+\alpha_B)$', 'color', colours_cb(5));
    plot(vs.param_1_vals(range), vs.alpha_B_vals(range).^2,...
        'c-','displayname', '$\alpha_B^2$', 'color', colours_cb(7));
    legend

elseif curve_to_plot==5
    % for bending energy in A
    plot(vs.param_1_vals(range), vs.E_all(5,range)./(vs.d_vals(range).^2*vs.kD),...
        'k-','displayname', '$E_\mathrm{Bend,A}$');
    plot(vs.param_1_vals(range), 4*pi*(1-cos(vs.phi_vals(range))),...
        '-','displayname', '$4 \pi (1-\cos{\phi})$', 'color', colours_cb(1));
    plot(vs.param_1_vals(range), vs.kappa_vals(range)./(vs.d_vals(range).^2*vs.kD),...
        'r-','displayname', '$\bar{\kappa}$', 'color', colours_cb(4));
    legend

elseif curve_to_plot==6
    % for bending energy in B
    plot(vs.param_1_vals(range), vs.E_all(6,range)./(vs.d_vals(range).^2*vs.kD),...
        'k-','displayname', '$E_\mathrm{Bend,B}$');
    plot(vs.param_1_vals(range), pi/2*sin(vs.phi_vals(range)).^4./cos(vs.phi_vals(range)).^2,...
        '-','displayname', '$\pi \sin^4{\phi}/2 \cos^2{\phi}$', 'color', colours_cb(1));
    plot(vs.param_1_vals(range), vs.sigma_vals(range),...
        'r-','displayname', '$\sigma$', 'color', colours_cb(4));
    rl1 = vs.R_vals.*sin(vs.phi_vals)./vs.lambda_vals;
    bes_funcs = (-besselk(0,rl1,1).^2+besselk(1,rl1,1).^2)./besselk(1,rl1,1).^2;
    plot(vs.param_1_vals(range), bes_funcs(range),...
        'g-','displayname', '$(K_1^2-K_0^2)/K_1^2$', 'color', colours_cb(5));
    plot(vs.param_1_vals(range), vs.alpha_B_vals(range),...
        'c-','displayname', '$\alpha_B$', 'color', colours_cb(7));
    legend
end

end