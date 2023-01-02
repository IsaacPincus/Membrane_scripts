% first get the initial unwrapped shape from the BSplines

% n = 68;
% s0 = 0:0.001:pi;
% knot_s0(1:3) = s0(1);
% knot_s0(4:n+2) = linspace(s0(1), s0(end), n-1);
% knot_s0(n+3:n+5) = s0(end);
% 
% psi_at_knots(1:3) = 0;
% psi_at_knots(4:n) = linspace(0,pi,n-3);
% 
% t = knot_s0;
% a = psi_at_knots;
% 
% figure();
% hold on
% fnplt(spmak(t,a))

sc0 = pi/2;
n = 68;
n_int = 1000;
L0 = 2*pi;
R0 = 2;

s0_3 = linspace(0,sc0,n);
% psi_3 = linspace(0,sc0,n);
psi_3 = s0_3/R0;
psi_dot0 = 1/R0;
psi_dot1 = 1/R0;

lam_s_3 = 2*ones(size(psi_3));

pp_psi_3 = csape(s0_3,[psi_dot0, psi_3, psi_dot1], 'clamped');

trapz(s0_3, psi_3);
int = integrate_pp(pp_psi_3);
int(end);

pp_lam_s = csape(s0_3, lam_s_3);
s_3 = integrate_pp(pp_lam_s);

% s0_1_f

s_3_fine = linspace(s_3(1), s_3(end), n_int);
s0_3_fine = linspace(s0_3(1), s0_3(end), n_int);
r_3 = cumtrapz(s_3_fine, cos(fnval(pp_psi_3, s0_3_fine)));
z_3 = cumtrapz(s_3_fine, sin(fnval(pp_psi_3, s0_3_fine)));

figure();
hold on
axis equal
plot(r_3,z_3)

%% get exterior points

s0_1 = linspace(sc0,L0,n);
psi_1 = s0_1/R0;
% s0_1 = linspace(0,sc0,n);
% psi_1 = linspace(0,sc0,n);
psi_dot0 = 1/R0;
psi_dot1 = 1/R0;

lam_s_1 = 2*ones(size(psi_1));

pp_psi_1 = csape(s0_1,[psi_dot0, psi_1, psi_dot1], 'clamped');

trapz(s0_1, psi_1);
int = integrate_pp(pp_psi_1);
int(end);

pp_lam_s = csape(s0_1, lam_s_1);
s_1 = integrate_pp(pp_lam_s);

% s0_1_f

s_1_fine = linspace(s_1(1), s_1(end), n_int);
s0_1_fine = linspace(s0_1(1), s0_1(end), n_int);
r_1 = r_3(end) + cumtrapz(s_1_fine, cos(fnval(pp_psi_1, s0_1_fine)));
z_1 = z_3(end) + cumtrapz(s_1_fine, sin(fnval(pp_psi_1, s0_1_fine)));

% figure();
% hold on
% axis equal
plot(r_1,z_1)

%% now final section

tc = s_3(end);
t_2 = linspace(tc,10,n);
t_3 = s_3;
psi_m_3 = psi_3;
% psi_3 = linspace(0,sc0,n);
psi_m_2 = zeros(size(t_2));
psi_dot0 = 1/R0;
psi_dot1 = 1/R0;

pp_psi_m_2 = csape(t_2,[psi_dot0, psi_m_2, psi_dot1], 'clamped');

t_2_fine = linspace(t_2(1), t_2(end), n_int);

r_2 = r_3(end) + cumtrapz(t_2_fine, cos(fnval(pp_psi_m_2, t_2_fine)));
z_2 = z_3(end) + cumtrapz(t_2_fine, sin(fnval(pp_psi_m_2, t_2_fine)));

plot(r_2,z_2)

%% get the energy function 



%% old stuff

% figure();
% hold on
% % fnplt(pp_psi_1)
% for ii=1:n-1
%     x = pp_psi_3.breaks(ii:ii+1);
%     y = polyval(pp_psi_3.coefs(ii,:),[0, x(2)-x(1)]);
%     plot(x, y)
% end
% plot(s0_3,psi_3,'ro','MarkerFaceColor','r')




