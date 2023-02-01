% implements Baulin 2021 microplastics model
% all SI units
clear variables

% sigma_vals=linspace(1e-4, 0.1, 20);
sigma_vals = logspace(-5,-3,20);
R = 0.015;
d = 2;
% sigma_vals = R^2/d^2;
alpha_i = 0;
% for zeta = [-0.005, -0.0025, -0.0005]
for zeta = [-0.001]
    for ii = 1:length(sigma_vals)
        sigma = sigma_vals(ii);
%         zeta = -0.01;
        
        R = sqrt(sigma*d^2);
    %     R = 0.2;
    %     theta = pi/3;
    %     S_i = d^2;
    %     S_A = 2*pi*R^2*(1-cos(theta));
    %     S_B = d^2 - pi*R^2*sin(theta)^2;
        
%         lambda = 0;
    %     sigma(ii) = R^2/d^2;
        
        % const = [zeta, alpha_i, S_i, S_A, S_B, lambda];
        % const = [zeta, alpha_i, S_i, S_A, S_B];
        
        % free([1.2,1], const)
        
        % alpha_A = linspace(-0.05,0.05, 100);
        % alpha_B = linspace(-0.05,0.05, 100);
        % 
        % [X,Y] = meshgrid(alpha_A, alpha_B);
        
        % surf(X, Y, free_surf(X, Y, lambda, const))
        % hold on
        
        options = optimoptions('fmincon','MaxFunEvals', 1e5,'MaxIter', 1e5,...
            'algorithm', 'sqp',...
            'OptimalityTolerance', 1e-9,'StepTolerance', 1e-9);
        
        % contour(X, Y, free_surf(X, Y, lambda, const), 50)
        % fimplicit(@(x,y) S_A./(1+x)+S_B./(1+y)-S_i/(1+alpha_i), [-0.05,0.05], 'r-')
        
        const = [zeta, alpha_i, d, R];
        % out = fminunc(@(y) free_theta(y, const), [0.005, 0.008, pi/3, 0.1])
        % out = fmincon(@(y) free(y, const),[0.005, 0.008, pi/3], [],[],[],[],...
        %     [-1,-1,0],[Inf, Inf, pi/2], ...
        %     @(y) area_con(y,const))
        alpha_A_init = 0;
        theta_init = deg2rad(20);
        const_constraint = [zeta, alpha_i, d, R, alpha_A_init, theta_init];
        alpha_B_init = fzero(@(y) area_con_theta_init(y,const_constraint), 0);
        [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
            fmincon(@(y) free_theta(y, const),[alpha_A_init, alpha_B_init, theta_init],...
            [],[],[],[],...
            [-0.1,-0.1,0],[0.1, 0.1, pi], ...
            @(y) area_con_theta(y,const), options);
        % % out = fminsearch(@(y) free(y, const), [1, 1])
        % 
        % % free(out,const)
        S_i = d^2;
        S_A = 2*pi*R^2*(1-cos(out(3)));
        S_B = d^2 - pi*R^2*sin(out(3))^2;
    %     (S_A/(1+out(1))+S_B/(1+out(2))-S_i/(1+alpha_i))
    %     zeta*S_A/(1+out(1)) + 1/2*(out(1)^2*S_A/(1+out(1))+out(2)^2*S_B/(1+out(2)))
        
        alpha_A(ii) = out(1);
        alpha_B(ii) = out(2);
        theta(ii) = out(3);

        lambda = lam_vals.eqnonlin;

        alpha_A_test(ii) = sqrt(1+2*(lambda+zeta))-1;
        alpha_B_test(ii) = sqrt(1+2*lambda)-1;

    end

    % figure();
    hold on
%     plot([theta,theta], [-1,1], ':')
    plot(sigma_vals, alpha_B, 'o', 'displayname', sprintf('$\\zeta=%g$', zeta))
%     plot(sigma_vals, alpha_A, 'x', 'displayname', sprintf('$\\zeta=%g$', zeta))
end


%% plot Baulin 2021 data

opts = detectImportOptions("Baulin_Fig1C.csv");
data = readmatrix('Baulin_Fig1C.csv', opts);
% opts = detectImportOptions("Baulin_Fig1D.csv");
% data = readmatrix('Baulin_Fig1D.csv', opts);
X = data(:,1:2:end);
Y = data(:,2:2:end);
% Y = data(:,2:2:end)*1.5902;

% axes1 = gca;
% axes1.XScale = 'log';

% plot(X, Y, 'displayName', 'Baulin 2021')
% legend

rad2deg(theta)

%% functions

function f = free(y, const)
    % function of free energy to minimise

    zeta = const(1);
    alpha_i = const(2);
    S_i = const(3);
    S_A = const(4);
    S_B = const(5);
    lambda = const(6);

    alpha_A = y(1);
    alpha_B = y(2);
    theta = y(3);
    S_A = 2*pi*R^2*(1-cos(theta));
    S_B = d^2 - pi*R^2*sin(theta)^2;
%     lambda = y(3);

    f = zeta*S_A./(1+alpha_A) ...
      + 1/2*(alpha_A.^2*S_A./(1+alpha_A)+alpha_B.^2*S_B./(1+alpha_B))...
      + lambda.*(S_A./(1+alpha_A)+S_B./(1+alpha_B)-S_i/(1+alpha_i));

end

function f = free_theta(y, const)
    % function of free energy to minimise

    zeta = const(1);
    alpha_i = const(2);
    d = const(3);
    R = const(4);

    alpha_A = y(1);
    alpha_B = y(2);
    theta = y(3);

    S_A = 2*pi*R^2*(1-cos(theta));
    S_B = d^2 - pi*R^2*sin(theta)^2;

    f = zeta*S_A./(1+alpha_A) ...
      + 1/2*(alpha_A.^2*S_A./(1+alpha_A)+alpha_B.^2*S_B./(1+alpha_B));

end

function f = free_theta_lambda(y, const)
    % function of free energy to minimise

    zeta = const(1);
    alpha_i = const(2);
    d = const(3);
    R = const(4);

    alpha_A = y(1);
    alpha_B = y(2);
    theta = y(3);
    lambda = y(4);

    S_i = d^2;
    S_A = 2*pi*R^2*(1-cos(theta));
    S_B = d^2 - pi*R^2*sin(theta)^2;

    f = zeta*S_A./(1+alpha_A) ...
      + 1/2*(alpha_A.^2*S_A./(1+alpha_A)+alpha_B.^2*S_B./(1+alpha_B))...
      + lambda.*(S_A./(1+alpha_A)+S_B./(1+alpha_B)-S_i/(1+alpha_i));

end

function f = free_surf(alpha_A, alpha_B, lambda, const)
    % function of free energy to minimise

    zeta = const(1);
    alpha_i = const(2);
    S_i = const(3);
    S_A = const(4);
    S_B = const(5);
%     lambda = const(6);

%     alpha_A = y(1);
%     alpha_B = y(2);
%     lambda = y(3);

    f = zeta*S_A./(1+alpha_A) ...
      + 1/2*(alpha_A.^2*S_A./(1+alpha_A)+alpha_B.^2*S_B./(1+alpha_B))...
      + lambda.*(S_A./(1+alpha_A)+S_B./(1+alpha_B)-S_i/(1+alpha_i));

end

function [c,ceq] = area_con(y, const)

    zeta = const(1);
    alpha_i = const(2);
    S_i = const(3);
    S_A = const(4);
    S_B = const(5);

    alpha_A = y(1);
    alpha_B = y(2);

    c(1) = 0;

    ceq(1) = S_A./(1+alpha_A)+S_B./(1+alpha_B)-S_i/(1+alpha_i);

end

function [c,ceq] = area_con_theta(y, const)

    zeta = const(1);
    alpha_i = const(2);
    d = const(3);
    R = const(4);

    alpha_A = y(1);
    alpha_B = y(2);
    theta = y(3);

    S_i = d^2;
    S_A = 2*pi*R^2*(1-cos(theta));
    S_B = d^2 - pi*R^2*sin(theta)^2;

    c(1) = 0;

    ceq(1) = S_A./(1+alpha_A)+S_B./(1+alpha_B)-S_i/(1+alpha_i);

end

function ceq = area_con_theta_init(y, const)

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