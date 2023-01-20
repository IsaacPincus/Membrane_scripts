N = 1e5;
d = 10;
R = 1;
epsilon = -1;
n0 = 1;
Sigma = 0.1;
kappa = 1;
lambda = sqrt(kappa/Sigma);

% get the shape of the free region
% phi_vals = linspace(0.001, pi/8, 10);
% h0_vals = linspace(-1.5,0.5,20);
phi_vals = deg2rad(60);
h0_vals = -R/2;

E = zeros(length(phi_vals), length(h0_vals));

figure();
hold on
% xlabel('$h_0$')
% ylabel('$E$')
axis equal

for ii = 1:length(phi_vals)
for jj = 1:length(h0_vals)

phi = phi_vals(ii);
h0 = h0_vals(jj);
r_phi = sin(phi)*R;
r = linspace(r_phi, d/2,N);
S_A = 2*pi*R^2*(1-cos(phi));

% [h,C,A] = free_shape_linear(r, R, d, phi, kappa, Sigma, -epsilon*n0);
[h,C,A,E(ii,jj),hderiv,lap_h] = free_shape_linear_fixed_h(r, r_phi, d, phi, kappa, Sigma, h0);
[C2,A2,E2(ii,jj)] = free_shape_linear_no_curve(r_phi, d, phi, kappa, Sigma, h0);

[r_nonlin, h_nonlin, out] = ...
    free_shape_nonlinear_fixed_h(R, d, phi, kappa, Sigma, h0, -1, -0.1);
out.y(1,end)
h_nonlin(end)

solution = deval(out, linspace(0,out.xe, 1000));

curves = solution(4,:).^2./(2*r_nonlin);

out.ye(7)*pi*kappa
kappa*pi*trapz(r_nonlin, curves)
kappa/2*2*pi*trapz(r, r.*lap_h.^2)

% 
% p_psi_1 = out.y(4,1)
% p_r_1 = out.y(5,1)

% E(ii,jj) = E(ii,jj) + epsilon*n0*S_A;

% hderiv = C(1)./r+C(3)/lambda*besseli(1,r/lambda)-C(4)/lambda*besselk(1,r/lambda);
% lap_h = C(3)/lambda^2*besselk(0,r/lambda)+C(4)/lambda^2*besselk(0,r/lambda);

% area_func = @(x) x.*sqrt(1+(C(1)./x+C(3)/lambda*besseli(1,x/lambda)-C(4)/lambda*besselk(1,x/lambda)).^2);
% 
% A_test = 2*pi*trapz(r, r.*sqrt(1+(hderiv).^2)) + d^2*(1-pi/4)
% A_test2 = 2*pi*integral(area_func,r_phi,d/2) + d^2*(1-pi/4)

% plot of shape and nanoparticle
l1 = plot(r,h, '-', 'DisplayName','linear');
plot(r_nonlin, h_nonlin, '-', 'DisplayName','nonlinear')
color = get(l1, 'color');
t = linspace(-pi/2,pi/2,1000);
x = cos(t)*R;
% y = sin(t)*R+(R*cos(phi)+h(1));
y = sin(t)*R+R*cos(phi)+h(1);
plot(x,y, '--', 'color', color, 'HandleVisibility','off')
% plot(x, (x-sin(phi))*tan(phi)+h(1))
annotation('textbox', [0.3,0.7,0.4,0.2], 'String',...
    [sprintf('$d = %0.2g$ \n', d),...
    sprintf('$R = %0.2g$ \n', R),...
    sprintf('$w = %0.2g$ \n', epsilon*n0),...
    sprintf('$\\kappa = %0.2g$ \n', kappa)], ...
    'interpreter', 'latex','FontSize',16,...
    'FitBoxToText','on','LineStyle','none', 'Color','k')
annotation('textbox', [0.55,0.7,0.4,0.2], 'String',...
    [sprintf('$\\Sigma = %0.2g$ \n', Sigma),...
    sprintf('$\\phi = %0.2g^{\\circ}$ \n', rad2deg(phi)),...
    sprintf('$h_0 = %0.2g$ \n', h0)], ...
    'interpreter', 'latex','FontSize',16,...
    'FitBoxToText','on','LineStyle','none', 'Color','k')

end
% plot(h0_vals, E(ii,:), 'displayname', sprintf('$\\phi=%0.2g$', phi));

end
legend('Location','southeast')

% E_min_const_h0 = min(E,[],1);
% E_min_const_phi = min(E,[],2);
% 
% figure();
% hold on
% xlabel('$h_0$')
% ylabel('$\phi$')
% zlabel('$E$')
% surf(h0_vals, phi_vals,  E)
% % plot(h0_vals, E);
% 
% figure();
% hold on
% xlabel('$\phi$')
% ylabel('$E$')
% plot(phi_vals, E_min_const_phi)
% 
% figure();
% hold on
% xlabel('$h_0$')
% ylabel('$E$')
% plot(h0_vals, E_min_const_h0)

% plot of derivative
% figure();
% hold on
% plot(r,hderiv)
% plot(r, tan(phi)*ones(size(r)))