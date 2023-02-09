clear variables
close all

N = 1e5;
d = 10;
R = 1;
epsilon = -1;
n0 = 1;
Sigma = 10;
kappa = 1;
lambda = sqrt(kappa/Sigma);

% get the shape of the free region
phi_vals = deg2rad(linspace(5, 160, 5));
% h0_vals = linspace(-1.5,0.5,20);
% phi_vals = deg2rad(160);
h_phi_vals = -R/2;
h_phi_init = -0.45;
psi_dot_init = 1;

E = zeros(length(phi_vals), length(h_phi_vals));

figure();
hold on
% xlabel('$h_0$')
% ylabel('$E$')
axis equal

for ii = 1:length(phi_vals)
    
    
    
    phi = phi_vals(ii);
    r_phi = sin(phi)*R;
    r = linspace(r_phi, d/2,N);
    S_A = 2*pi*R^2*(1-cos(phi));
    if phi<pi/2
        psi_dot_init = -1;
    else
        psi_dot_init = 1;
    end
    
    % [h,C,A] = free_shape_linear(r, R, d, phi, kappa, Sigma, -epsilon*n0);
    % [h,C,A,E(ii),hderiv,lap_h] = free_shape_linear_fixed_h(r, r_phi, d, phi, kappa, Sigma, h_phi_init);
    % [C2,A2,E2(ii,jj)] = free_shape_linear_no_curve(r_phi, d, phi, kappa, Sigma, h0);
    
    constants = [r_phi, d, phi, kappa, Sigma, N];
    options = optimoptions('fmincon', 'algorithm','sqp', 'MaxFunctionEvaluations',1000);
    h_phi = fmincon(@(inp) min_h_phi(inp, constants), h_phi_init,...
        [],[],[],[],[],[],[],options);
    [h,C,A,E(ii),hderiv,lap_h] = free_shape_linear_fixed_h(r, r_phi, d, phi, kappa, Sigma, h_phi);
    
    %%%%% normal full solution
    % [r_nonlin, h_nonlin, out] = ...
    %     free_shape_nonlinear_fixed_h(R, d, phi, kappa, Sigma, h0, -1, -0e.1);
    out = free_shape_nonlinear_free_h(R, d, phi, kappa, Sigma, 1);
    % out.y(1,end)
    % h_nonlin(end)
    p_psi_init_final = out.y(4,1);
    psi_dot_init_final(ii) = out.y(4,1)/(2*r_phi)-sin(out.y(1,1))/r_phi

    h_phi = 0;
    
    solution = deval(out, linspace(0,out.xe, 1000));
    r_nonlin = solution(2,:);
    h_nonlin = solution(3,:)-R*sin(phi)+h_phi;
    E_nonlin(ii) = solution(7,end)*kappa*pi;
    E_lin(ii) = kappa/2*2*pi*trapz(r, r.*lap_h.^2);
    S_B_lin(ii) = A;
    A_nonlin = solution(8,end)*2*pi;
    if out.x(end)<d/2
        r_nonlin = [solution(2,1:end-1),d/2];
        h_nonlin = [solution(3,1:end-1),solution(3,end)]-R*sin(phi)+h_phi;
        A_nonlin = A_nonlin+pi*((d/2)^2-solution(2,end)^2);
    elseif out.x(end)>d/2
        r_nonlin(r_nonlin>d/2) = d/2;
        h_nonlin(r_nonlin>d/2) = solution(3,end)-R*sin(phi)+h_phi;
        A_nonlin = A_nonlin+pi*((d/2)^2-solution(2,end)^2);
    end
    S_B_nonlin(ii) = A_nonlin + d^2*(1-pi/4);
    
    curves = solution(4,:).^2./(2*r_nonlin);
    %%%%%%
    
    %%%%%%% single solutions
    % const(1) = lambda;
    % const(2) = d;
    % 
    % % psi_dot_init_final = -1;
    % psi_dot_init_final = -1;
    % p_psi_init_final = 2*r_phi*(psi_dot_init_final+sin(phi)/r_phi);
    % p_r_init_final = -(p_psi_init_final^2/(4*r_phi)-p_psi_init_final*sin(phi)/r_phi...
    %     -2*r_phi/lambda^2)/cos(phi);
    % 
    % init_vals = [phi, r_phi, R*sin(phi), p_psi_init_final, p_r_init_final, 0,0];
    % event_func = @(s,y) myEventFcn(s,y,const);
    % %     options = odeset('Events', @myEventFcn);
    % options = odeset('Events', event_func);
    % out = ode45(@(s, y) hamilton(s, y,const),linspace(0,10,1000),init_vals,...
    %     options);
    % 
    % %     solution = deval(out, linspace(0,out.xe, 1000));
    % solution = deval(out, linspace(0,out.x(end), 1000));
    % r_nonlin = solution(2,:);
    % h_nonlin = solution(3,:)-R*sin(phi)+h_phi;
    % 
    % psi_end = out.y(1,end);
    % diffs = psi_end^2
    %%%%%%
    
    % out.y(7,end)*pi*kappa
    % kappa*pi*trapz(r_nonlin, curves)
    % kappa/2*2*pi*trapz(r, r.*lap_h.^2)
    
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
%     l1 = plot(r,h, '-', 'DisplayName','linear');
    l1 = plot(r_nonlin, h_nonlin, '-', 'DisplayName','nonlinear');
    color = get(l1, 'color');
    t = linspace(-pi/2,pi/2,1000);
    x = cos(t)*R;
%     y = sin(t)*R+(R*cos(phi)+h(1));
    % y = sin(t)*R+R*cos(phi)+h_phi;
    y = sin(t)*R+R*cos(phi)+h_nonlin(1);
    plot(x,y, '--', 'color', color, 'HandleVisibility','off')
%     plot(x,y, '--','HandleVisibility','off')
    xlabel('$r$')
    ylabel('$h$')
    % plot(x, (x-sin(phi))*tan(phi)+h(1))
    % annotation('textbox', [0.3,0.7,0.4,0.2], 'String',...
    %     [sprintf('$d = %0.2g$ \n', d),...
    %     sprintf('$R = %0.2g$ \n', R),...
    %     sprintf('$w = %0.2g$ \n', epsilon*n0),...
    %     sprintf('$\\kappa = %0.2g$ \n', kappa)], ...
    %     'interpreter', 'latex','FontSize',16,...
    %     'FitBoxToText','on','LineStyle','none', 'Color','k')
    % annotation('textbox', [0.55,0.7,0.4,0.2], 'String',...
    %     [sprintf('$\\Sigma = %0.2g$ \n', Sigma),...
    %     sprintf('$\\phi = %0.2g^{\\circ}$ \n', rad2deg(phi)),...
    %     sprintf('$h_0 = %0.2g$ \n', h_phi_init)], ...
    %     'interpreter', 'latex','FontSize',16,...
    %     'FitBoxToText','on','LineStyle','none', 'Color','k')

end
% plot(h0_vals, E(ii,:), 'displayname', sprintf('$\\phi=%0.2g$', phi));

% legend('Location','southeast')

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

function E = min_h_phi(h_phi_guess, constants)

    r_phi = constants(1);
    d = constants(2);
    phi = constants(3);
    kappa = constants(4);
    Sigma = constants(5);
    N = constants(6);

    r = linspace(r_phi, d/2,N);
    [~,~,~,E,~,~] = free_shape_linear_fixed_h(r, r_phi, d, phi, kappa, Sigma, h_phi_guess);
end

function derivs = hamilton(s, y, const)
    psi = y(1);
    r = y(2);
    h = y(3);
    p_psi = y(4);
    p_r = y(5);
    p_h = y(6);
    E_bend = y(7);

    lambda = const(1);
    d = const(2);

    f(1) = p_psi/(2*r)-sin(psi)/r;
    f(2) = cos(psi);
    f(3) = sin(psi);
    f(4) = cos(psi)*(p_psi/r-p_h)...
        +sin(psi)*p_r;
    f(5) = p_psi/r*(p_psi/(4*r)-sin(psi)/r)...
        +2/lambda^2;
    f(6) = 0;
    f(7) = p_psi^2/(2*r);

    derivs = f';
end

function [value,isterminal,direction] = myEventFcn(s,y,const)
    psi = y(1);
    r = y(2);
    h = y(3);
    p_psi = y(4);
    p_r = y(5);
    p_h = y(6);

    lambda = const(1);
    d = const(2);

%     value = y(2)-d/2;
%     value = y(1);
%     isterminal = 1;
%     direction = 0;

    value = [y(1);y(1)-pi];
    isterminal = [1;1];
    direction = [0;0];

end