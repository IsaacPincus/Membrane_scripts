function fig_handle = membrane_shapes(xaxis_name, xscale, yaxis_name, yscale,...
    vs, patch_lims, plot_curves)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

fig_handle = figure('Position',[400,100,800,600]);
hold on
axis equal
xlabel('$r$ ($\mu$m)')
ylabel('$h$ ($\mu$m)')
% ylabel('$\int r \nabla^2 h dr$')
% axes1 = gca;
% axes1.YScale = 'log';
N = 1e4;
counter = 0;
for ii = plot_curves
    counter = counter+1;

    alpha_A = vs.alpha_A_vals(ii);
    alpha_B = vs.alpha_B_vals(ii);
    phi = vs.phi_vals(ii);
    h_phi = vs.h_phi_vals(ii);
    
    % get the shape of the free region
    Sigma = vs.Sigma_vals(ii);
    lambda = vs.lambda_vals(ii);
    d = vs.d_vals(ii);
    kappa = vs.kappa_vals(ii);
    epsilon = vs.epsilon;
    kD = vs.kD;
    R = vs.R_vals(ii);
    n0 = vs.n0;
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
%     h1 = plot(r, h, 'displayname', sprintf('$\\kappa = %0.2e$', ...
%         kappa), 'Color', colours_cb(counter));
    h1 = plot(r, h, 'displayname', sprintf('$\\kappa = %0.2e$', ...
        kappa), 'Color', 'k');
%     h2 = plot(r, h_lim, '--', 'displayname', sprintf('$R = %0.2e$', ...
%         R), 'Color', colours_cb(counter), 'HandleVisibility','off');
%     yyaxis left
%     h1 = plot(r, r.*lap_h.^2,'-', 'displayname', sprintf('$R = %0.2e$', ...
%         R), 'Color', colours_cb(counter));
%     yyaxis right
%     h2 = plot(r, cumtrapz(r,r.*lap_h.^2),'-',  'displayname', sprintf('$R = %0.2e$', ...
%         R), 'Color', colours_cb(counter));
    result(counter) = trapz(r,r.*lap_h.^2);
    rl = r_phi/lambda;
    result2(counter) = tan(phi)^2/besselk(1,rl,1)^2*rl^2*1/2*(-besselk(0,rl,1)^2+besselk(1,rl,1)^2);
    colour = h1.Color;
    t = linspace(-pi/2,-pi/2+phi,1000);
%     t = linspace(-pi/2,pi/2,1000);
    x = cos(t)*R;
    % y = sin(t)*R+(R*cos(phi)+h(1));
    y = sin(t)*R+R*cos(phi)+h(1);
    plot(x,y, ':', 'HandleVisibility','off', 'color', colour)
end

legend

end