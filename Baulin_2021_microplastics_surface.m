% implements Baulin 2021 microplastics model
% all SI units
clear variables

d = 2;
R = 0.01;
alpha_i = 0;
zeta = -0.001;
sigma = R^2/d^2;

N = 100;
alpha_A_vals = linspace(-0.1,0.1,N);
theta_vals = linspace(0,pi,N);
Z = zeros(N,N);

const = [zeta, alpha_i, d, R];
for aa=1:length(alpha_A_vals)
    for tt = 1:length(theta_vals)
        alpha_A = alpha_A_vals(aa);
        theta = theta_vals(tt);
        Z(tt,aa) = free_theta(alpha_A, theta, const);
    end
end

figure();
surf(alpha_A_vals, theta_vals, Z);

[M,I] = min(Z, [], "all");
ai = ceil(I/N);
ti = mod(I,N);
alpha_A_opt = alpha_A_vals(ai)
theta_opt = rad2deg(theta_vals(ti))


%% functions

function f = free_theta(alpha_A, theta, const)
    % function of free energy to minimise

    zeta = const(1);
    alpha_i = const(2);
    d = const(3);
    R = const(4);

    const_constraint = [zeta, alpha_i, d, R, alpha_A, theta];
    alpha_B = fzero(@(y) area_con_theta(y,const_constraint), 0);

    S_A = 2*pi*R^2*(1-cos(theta));
    S_B = d^2 - pi*R^2*sin(theta).^2;

    f = zeta*S_A./(1+alpha_A) ...
      + 1/2*(alpha_A.^2*S_A./(1+alpha_A)+alpha_B.^2*S_B./(1+alpha_B));

end

function ceq = area_con_theta(y, const)

    zeta = const(1);
    alpha_i = const(2);
    d = const(3);
    R = const(4);
    alpha_A = const(5);
    theta = const(6);

    alpha_B = y;

    S_i = d^2;
    S_A = 2*pi*R^2*(1-cos(theta));
    S_B = d^2 - pi*R^2*sin(theta).^2;

    ceq(1) = S_A./(1+alpha_A)+S_B./(1+alpha_B)-S_i/(1+alpha_i);

end