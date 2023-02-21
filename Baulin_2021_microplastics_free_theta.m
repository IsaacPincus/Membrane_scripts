% implements Baulin 2021 microplastics model
% all SI units
clear variables

% constants
R = 0.5;
sigma = 3.5e-5;
d = sqrt(R^2/sigma);    % um
% d = 20;    % um
phi = pi/12;
kD  = 300/10^12*1e9;    % picoJ/um^2
zeta = -2.7e-4;            % dimensionless
alpha_i = 0.01;

% d_vals = linspace(1.8,1.7,4);
% for d = d_vals
% d = 3e2;
% R = 1;
% alpha_i = 0.002;
% zeta = -0.1*alpha_i^2/2;
% sigma = R^2/d^2;
% theta_vals = flip(deg2rad(linspace(0.01,5,50)));
theta_vals = deg2rad(linspace(0.01,30,50));
% theta_vals = deg2rad(30);
% theta_vals = deg2rad(4.898163265306122);
for ii=1:length(theta_vals)

    theta = theta_vals(ii);

%     options = optimoptions('fmincon','MaxFunEvals', 1e5, 'MaxIter', 1e3, 'algorithm', 'interior-point',...
%         'OptimalityTolerance', 1e-6, 'ConstraintTolerance', 1e-6, 'StepTolerance', 1e-6,...
%         'SpecifyObjectiveGradient', true, 'CheckGradients', true);

%     options = optimoptions('fmincon','MaxFunEvals', 1e5, 'MaxIter', 1e3, 'algorithm', 'sqp',...
%         'OptimalityTolerance', 1e-6, 'ConstraintTolerance', 1e-6, 'StepTolerance', 1e-6,...
%         'SpecifyObjectiveGradient', true, 'CheckGradients', true,...
%         'FiniteDifferenceType', 'central', 'FiniteDifferenceStepSize', 1e-10,...
%         'SpecifyConstraintGradient', true);

%     options = optimoptions('fmincon','MaxFunEvals', 1e5, 'MaxIter', 1e3, 'algorithm', 'sqp',...
%         'OptimalityTolerance', 1e-12, 'ConstraintTolerance', 1e-12, 'StepTolerance', 1e-12,...
%         'SpecifyObjectiveGradient', true,...
%         'SpecifyConstraintGradient', true);

%     options = optimoptions('fmincon','MaxFunEvals', 1e5, 'MaxIter', 1e3, 'algorithm', 'sqp',...
%             'OptimalityTolerance', 1e-12, 'ConstraintTolerance', 1e-12, 'StepTolerance', 1e-12,...
%             'Diagnostics', 'on',...
%             'Display', 'iter-detailed');
% 
%     const = [zeta, alpha_i, d, R, theta];
%     [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
%         fmincon(@(y) free_theta(y, const),[(alpha_i+zeta)/(1-zeta), alpha_i], [],[],[],[],...
%         [-0.1,-0.1],[0.1, 0.1], ...
%         @(y) area_con_theta(y,const), options);

    clear options
    options.verify_level = 1;
    options.printfile = 'snsolve_print.txt';
    options.scale_option = 2;
    options.major_feasibility_tolerance = 1e-12;
    options.major_optimality_tolerance = 1e-12;
%     options.function_precision = 1e-8;
%     options.major_optimality_tolerance = 1e-4;
    const = [zeta, alpha_i, d, R, theta];
    [out,fval,exitflag,output,lam_vals,states] = ...
        snsolve(@(y) free_theta(y, const),[(alpha_i+zeta)/(1-zeta), alpha_i], [],[],[],[],...
        [-0.1,-0.1]',[0.1, 0.1]', ...
        @(y) area_con_theta(y,const), options);
    
    alpha_A = out(1);
    alpha_B = out(2);

    S_A = 2*pi*R^2*(1-cos(theta));
    S_B = d^2 - pi*R^2*sin(theta)^2;

    E_vals(ii) = fval;
    alpha_A_vals(ii) = out(1);
    alpha_B_vals(ii) = out(2);

end

% initial_E = alpha_i^2*d^2/(1+alpha_i)/2;
initial_E = 0;

% figure();
hold on
plot(rad2deg(theta_vals), (E_vals-initial_E)*kD);
% ylim([min(E_vals), max(E_vals)])
xlabel('$\phi$')
ylabel('$E_\mathrm{total}$')

% alpha_i^2/2
% zeta

% end


%% functions

function [f, g] = free_theta(y, const)
    % function of free energy to minimise

    zeta = const(1);
    alpha_i = const(2);
    d = const(3);
    R = const(4);
    theta = const(5);

    alpha_A = y(1);
    alpha_B = y(2);

    S_A = 2*pi*R^2*(1-cos(theta));
    S_B = d^2 - pi*R^2*sin(theta)^2;

    g = [-zeta*S_A./(1+alpha_A)^2+alpha_A*S_A/(1+alpha_A)-1/2*alpha_A^2*S_A/(1+alpha_A)^2;
         alpha_B*S_B/(1+alpha_B)-1/2*alpha_B^2*S_B/(1+alpha_B)^2];

%         parts(1) = -kD*1/2*d^2*alpha_B^2/(1+alpha_B)^2;
%         parts(2) = kD*d^2*alpha_B/(1+alpha_B);
%         parts(3) = kD/2*pi*r_phi^2*alpha_B^2/(1+alpha_B)^2;
%         parts(4) = -kD*pi*r_phi^2*alpha_B/(1+alpha_B);

    f = zeta*S_A./(1+alpha_A) ...
      + 1/2*(alpha_A.^2*S_A./(1+alpha_A)+alpha_B.^2*S_B./(1+alpha_B))...
      - alpha_i^2*d^2/(1+alpha_i)/2;

end

function [c,ceq, gc, gceq] = area_con_theta(y, const)

    zeta = const(1);
    alpha_i = const(2);
    d = const(3);
    R = const(4);
    theta = const(5);

    alpha_A = y(1);
    alpha_B = y(2);

    S_i = d^2;
    S_A = 2*pi*R^2*(1-cos(theta));
    S_B = d^2 - pi*R^2*sin(theta)^2;

    c = [];

    ceq = S_A./(1+alpha_A)+S_B./(1+alpha_B)-S_i/(1+alpha_i);

    gc = [];

    gceq = [-S_A/(1+alpha_A)^2  -S_B/(1+alpha_B)^2]



end