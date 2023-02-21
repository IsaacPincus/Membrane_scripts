clear variables

MyTaskID = 0;
NumberOfTasks = 1;

% check that the environment variables have been read in correctly
if ~(exist('MyTaskID', 'var')&&exist('NumberOfTasks', 'var'))
    error('Environment variables not set correctly')
end

% taskIDs count from zero, alter this here
MyTaskID = MyTaskID + 1;

% generate list of independent variables to run, which should be in the
% order [epsilon, n0, d, R, kD, kappa, alpha_i] for each row
R_vals = 0.5;                 % um
sigma_vals = logspace(-4,-1,30);                % surface fraction
%d_vals = sqrt(R_vals.^2./sigma_vals);           % um
kD_vals = 300/10^12*1e9;                        % picoJ/um^2
kD_base = 300/10^12*1e9;                        % picoJ/um^2
zeta_vals = logspace(-6,-1,6);                  % dimensionless
epsilon_vals = -zeta_vals*kD_base;              % picoJ/um^2
n0_vals = 1;                                    % fraction
kappa_vals = 1e-19*1e12;         % picoJ
% kappa_vals = logspace(-21,-15, 60)*1e12;
alpha_i_vals = 0;                            % fraction
% other constants
N = 3e5;

ii = 0;
for rr = 1:length(R_vals)
    for ss = 1:length(sigma_vals)
        for kk = 1:length(kD_vals)
            for ee = 1:length(epsilon_vals)
                for nn = 1:length(n0_vals)
                    for pp = 1:length(kappa_vals)
                        for aa = 1:length(alpha_i_vals)
                            ii = ii+1;
                            parameter_set(ii, 1:7) = [epsilon_vals(ee),...
                                n0_vals(nn),...
                                sqrt(R_vals(rr)^2/sigma_vals(ss)),...
                                R_vals(rr),...
                                kD_vals(kk),...
                                kappa_vals(pp),...
                                alpha_i_vals(aa)];
                        end
                    end
                end
            end
        end
    end
end


%% split up job between task IDs
size_parameter_set = size(parameter_set);
total_parameter_sets = size_parameter_set(1);
min_sets_per_job = floor(total_parameter_sets/NumberOfTasks);
extra_sets = mod(total_parameter_sets, NumberOfTasks);
if MyTaskID <= extra_sets
    my_set_min = (min_sets_per_job+1)*(MyTaskID-1)+1
    my_set_max = (min_sets_per_job+1)*(MyTaskID)
else
    my_set_min = (min_sets_per_job+1)*(extra_sets) ...
        + min_sets_per_job*(MyTaskID-extra_sets-1) + 1
    my_set_max = (min_sets_per_job+1)*(extra_sets) ...
        + min_sets_per_job*(MyTaskID-extra_sets)
end
parameter_set_self = parameter_set(my_set_min:my_set_max, :);
my_set_size = size(parameter_set_self);
my_set_length = my_set_size(1);

save_name = sprintf('data/task%i_results.mat', MyTaskID-1);

%% run batch job for this task set
for ii = 1:my_set_length
    epsilon = parameter_set_self(ii, 1);
    n0 = parameter_set_self(ii, 2);
    d = parameter_set_self(ii, 3);
    R = parameter_set_self(ii, 4);
    kD = parameter_set_self(ii, 5);
    kappa = parameter_set_self(ii, 6);
    alpha_i = parameter_set_self(ii, 7);
    zeta = epsilon*n0/kD;
    sigma = sqrt(R^2/d^2);

    % solve microplastics for initial stretch of A and B
    options = optimset('MaxFunEvals', 1e5, 'MaxIter', 1e4);
    const = [zeta, alpha_i, d, R];
    [out,fval,~,~,lam_vals,grad,hessian] = ...
        fmincon(@(y) free_phi(y, const),[0.2, 0.2, pi/3], [],[],[],[],...
        [-1,-1,0],[Inf, Inf, pi/2], ...
        @(y) lipid_con_phi(y,const), options);
    alpha_A_init = out(1);
    alpha_B_init = out(2);
    phi_init = out(3);
    S_A_init = 2*pi*R^2*(1-cos(phi_init));
    S_B_init = d^2 - pi*R^2*sin(phi_init)^2;

    E_adhesion_init = epsilon*n0*S_A_init ./(1+alpha_A_init );
    E_stretch_A_init  = kD/2*(alpha_A_init .^2*S_A_init ./(1+alpha_A_init ));
    E_stretch_B_init  = kD/2*(alpha_B_init .^2*S_B_init ./(1+alpha_B_init ));
    E_total_init = E_adhesion_init+E_stretch_A_init+E_stretch_B_init;

    constants = [epsilon, n0, d, R, kD, kappa, alpha_i];
    inputs = [alpha_A_init, alpha_B_init, phi_init];
   
    [out,fval2,exitflag,output,solutions] = ...
        get_linear_minimum(constants, inputs);

    solution_set{ii} = solutions;
    
    % get energies etc
    alpha_A = out(1);
    alpha_B = out(2);
    phi = out(3);

    Sigma = kD*alpha_B;
    r_phi = sin(phi)*R;
    lambda = sqrt(kappa/Sigma);
    lamt = lambda*sqrt(alpha_B);

    [~, delA, Ebend] =...
        free_shape_linear_free_h(R, phi, kappa, Sigma);

    S_A = 2*pi*R^2*(1-cos(phi));
    S_B = delA+(d^2-pi*r_phi^2);

    E_adhesion = epsilon*n0*S_A./(1+alpha_A);
    E_stretch_A = kD/2*(alpha_A.^2*S_A./(1+alpha_A));
    E_stretch_B = kD/2*(alpha_B.^2*S_B./(1+alpha_B));
    E_bend_B = Ebend;
    E_bend_A = 4*pi*kappa*(1-cos(phi));

    % stretching, adhesion and bending energy
    E = E_adhesion + E_stretch_A + E_stretch_B + E_bend_A + E_bend_B;
    
    E_all(1,ii) = E;
    E_all(2,ii) = E_adhesion;
    E_all(3,ii) = E_stretch_A;
    E_all(4,ii) = E_stretch_B;
    E_all(5,ii) = E_bend_A;
    E_all(6,ii) = E_bend_B;
    
    alpha_A_vals(ii) = alpha_A;
    alpha_B_vals(ii) = alpha_B;
    phi_vals(ii) = phi;
    S_A_vals(ii) = S_A;
    S_B_vals(ii) = S_B;
    Sigma_vals(ii) = Sigma;
    
end

save(save_name)

%% extra functions
function f = stretch_bend_min(y, const)
    % function of free energy to minimise

    epsilon = const(1);
    n0 = const(2);
    d = const(3);
    R = const(4);
    kD = const(5);
    kappa = const(6);
    alpha_i = const(7);
    N = const(8);
    
    alpha_A = y(1);
    alpha_B = y(2);
    phi = y(3);
    h_phi = y(4);

    % problems when phi = 0, so change slightly
    if phi==0
        phi = phi+eps(1);
    end
    
    Sigma = kD*alpha_B;
%     lambda = sqrt(kappa/Sigma);
    r_phi = sin(phi)*R;
%     r = linspace(r_phi, d/2,N);
    
    [~,S_B, ~,E_bend] = free_shape_linear_no_curve(...
        r_phi, d, phi, kappa, Sigma, h_phi);
    
    S_A = 2*pi*R^2*(1-cos(phi));
    
    % stretching, adhesion and bending energy
    f = epsilon*n0*S_A./(1+alpha_A) ...
      + kD/2*(alpha_A.^2*S_A./(1+alpha_A)+alpha_B.^2*S_B./(1+alpha_B))...
      + E_bend + 4*pi*kappa*(1-cos(phi));

    [~, warnID] = lastwarn;
    if strcmp('MATLAB:integral:MaxIntervalCountReached', warnID)
        sprintf('%g,', y)
        warning('MATLAB:test_warning','Set new warning')
    end

end

function [c,ceq] = lipid_con_bend(y, const)

    epsilon = const(1);
    n0 = const(2);
    d = const(3);
    R = const(4);
    kD = const(5);
    kappa = const(6);
    alpha_i = const(7);
    N = const(8);
    
    alpha_A = y(1);
    alpha_B = y(2);
    phi = y(3);
    h_phi = y(4);

    % problems when phi = 0, so change slightly
    if phi==0
        phi = phi+eps(1);
    end
    
    Sigma = kD*alpha_B;
    r_phi = sin(phi)*R;
    
%     r = linspace(r_phi, d/2,N);
    [~,S_B, ~,~] = free_shape_linear_no_curve(...
        r_phi, d, phi, kappa, Sigma, h_phi);
    
    S_i = d^2;
    S_A = 2*pi*R^2*(1-cos(phi));

    c(1) = 0;

    ceq(1) = S_A./(1+alpha_A)+S_B./(1+alpha_B)-S_i/(1+alpha_i);

    [~, warnID] = lastwarn;
    if strcmp('MATLAB:integral:MaxIntervalCountReached', warnID)
        sprintf('%g,', y)
        warning('MATLAB:test_warning','Set new warning')
    end

end

function f = free_phi(y, const)
    % function of free energy to minimise

    zeta = const(1);
    alpha_i = const(2);
    d = const(3);
    R = const(4);

    alpha_A = y(1);
    alpha_B = y(2);
    phi = y(3);

    S_A = 2*pi*R^2*(1-cos(phi));
    S_B = d^2 - pi*R^2*sin(phi)^2;

    f = zeta*S_A./(1+alpha_A) ...
      + 1/2*(alpha_A.^2*S_A./(1+alpha_A)+alpha_B.^2*S_B./(1+alpha_B));

end

function [c,ceq] = lipid_con_phi(y, const)

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

function [out,fval2,exitflag,output,solutions] = ...
    get_linear_minimum(constants, inputs)

    epsilon = constants(1);
    n0 = constants(2);
    d = constants(3);
    R = constants(4);
    kD = constants(5);
    kappa = constants(6);
    alpha_i = constants(7);
    
    alpha_A_init = inputs(1);
    alpha_B_init = inputs(2);
    phi_init = inputs(3);

    alpha_A = alpha_A_init;
    alpha_B = alpha_B_init;
    phi = phi_init;
    
    Sigma = kD*alpha_B;
    r_phi = sin(phi)*R;
    lambda = sqrt(kappa/Sigma);
    lamt = lambda*sqrt(alpha_B);

    S_A = 2*pi*R^2*(1-cos(phi));
    S_B = (d^2-pi*r_phi^2);
    S_i = d^2;

%     options = optimoptions('fmincon','MaxFunEvals', 1e5, 'MaxIter', 1e3, 'algorithm', 'sqp',...
%         'OptimalityTolerance', 1e-12, 'ConstraintTolerance', 1e-12, 'StepTolerance', 1e-12,...
%         'SpecifyObjectiveGradient', true, 'CheckGradients', true);

%     options = optimoptions('fmincon','MaxFunEvals', 1e5, 'MaxIter', 1e3, 'algorithm', 'sqp',...
%         'OptimalityTolerance', 1e-12, 'ConstraintTolerance', 1e-12, 'StepTolerance', 1e-12,...
%         'SpecifyObjectiveGradient', true);

%     options = optimoptions('fmincon','MaxFunEvals', 1e5, 'MaxIter', 1e3, 'algorithm', 'sqp',...
%         'OptimalityTolerance', 1e-14, 'ConstraintTolerance', 1e-14, 'StepTolerance', 1e-14,...
%         'SpecifyConstraintGradient', true, 'SpecifyObjectiveGradient', true,...
%         'Diagnostics', 'on', 'Display', 'iter-detailed');

    options = optimoptions('fmincon','MaxFunEvals', 1e5, 'MaxIter', 1e3, 'algorithm', 'sqp',...
        'OptimalityTolerance', 1e-12, 'ConstraintTolerance', 1e-12, 'StepTolerance', 1e-12,...
        'SpecifyObjectiveGradient', false);

    problem = createOptimProblem("fmincon",...
        'x0', [alpha_A_init, alpha_B_init, phi_init],...
        'objective', @objective,...
        'lb', [-0.1,-0.1,0],...
        'ub', [0.1, 0.1, pi/2],...
        'nonlcon', @constraint,...
        'options', options);
%     alpha_A_vals_inp = [-0.1, -alpha_i, -1e-3, 0, 1e-3, alpha_i, 0.1];
%     phi_vals_inp = deg2rad([0, 5, 20, 45, 60, 85, 90]);
    alpha_A_vals_inp = [0, 1e-3, alpha_i, 0.1];
    phi_vals_inp = deg2rad([0,45,90]);
    ptmatrix_this_run =...
            zeros(length(alpha_A_vals_inp),length(phi_vals_inp), 3);
    tic
    for aa = 1:length(alpha_A_vals_inp)
        for pp = 1:length(phi_vals_inp)
            phi = phi_vals_inp(pp);
            alpha_B_vals_inp(pp,aa) =...
                fzero(@lipid_con_bend_test_A, rand()*0.2-0.1);
            ptmatrix_this_run(aa,pp,1:3) = ...
                [alpha_A_vals_inp(aa), alpha_B_vals_inp(pp,aa),...
                 phi_vals_inp(pp)];
        end
    end
    toc
    ptmatrix = reshape(ptmatrix_this_run, [numel(ptmatrix_this_run)/3, 3]);
    tpoints = CustomStartPointSet(ptmatrix);
    rs = RandomStartPointSet('NumStartPoints',10);
    gs = MultiStart("FunctionTolerance",1e-3, "XTolerance", 1e-3);
    tic
    [out,fval2,exitflag,output,solutions] = run(gs,problem, {tpoints,rs});
    toc
% 
%     clear options
%     options.verify_level = 1;
%     options.printfile = 'snsolve_print.txt';
%     options.scale_option = 2;
%     options.major_feasibility_tolerance = 1e-12;
%     options.major_optimality_tolerance = 1e-12;
% %     options.function_precision = 1e-8;
% %     options.major_optimality_tolerance = 1e-4;
%     [out,fval,exitflag,output,lam_vals,states] = ...
%         snsolve(@objective,[alpha_A_init, alpha_B_init], [],[],[],[],...
%         [-0.1,-0.1]',[0.1, 0.1]', ...
%         @constraint, options);

    function [f,g] = objective(y_obj)

        alpha_A = y_obj(1);
        alpha_B = y_obj(2);
        phi = y_obj(3);
        
        Sigma = kD*alpha_B;
        r_phi = sin(phi)*R;
        lambda = sqrt(kappa/Sigma);
        lamt = lambda*sqrt(alpha_B);

        [~, delA, Ebend] =...
            free_shape_linear_free_h(R, phi, kappa, Sigma);
    
        S_A = 2*pi*R^2*(1-cos(phi));
        S_B = delA+(d^2-pi*r_phi^2);

        alpha_A_gradient = (kD/2*S_A*alpha_A*(alpha_A+2)-epsilon*n0*S_A)...
            /(alpha_A+1)^2;

%         alpha_A_gradient = -epsilon*n0*S_A./(1+alpha_A)^2+kD*alpha_A*S_A/(1+alpha_A)...
%             -kD*1/2*alpha_A^2*S_A/(1+alpha_A)^2

%         parts(1) = -kD*1/2*d^2*alpha_B^2/(1+alpha_B)^2;
%         parts(2) = kD*d^2*alpha_B/(1+alpha_B);
%         parts(3) = kD/2*pi*r_phi^2*alpha_B^2/(1+alpha_B)^2;
%         parts(4) = -kD*pi*r_phi^2*alpha_B/(1+alpha_B);
%         parts(5) = -kD*pi*r_phi^2*alpha_B^2*tan(phi)^2/(1+alpha_B)^2/4;
%         parts(6) = +kD*3/4*r_phi^2*alpha_B*tan(phi)^2/(1+alpha_B)^2;
%         parts(7) = -pi*r_phi^2*kappa*tan(phi)^2/(2*lamt^2);
%         parts(8) = +(pi*besselk(0,r_phi/lambda,1)* ...
%         (-4*kD*r_phi*sqrt(alpha_B)*lamt^4 ...
%          -2*kD*r_phi*alpha_B^(3/2)*lamt^4 ...
%          -2*r_phi^3*sqrt(alpha_B)*kappa ...
%          -4*r_phi^3*alpha_B^(3/2)*kappa ...
%          -2*r_phi^3*alpha_B^(5/2)*kappa ...
%          +kD*r_phi^3 *alpha_B^(3/2)*lamt^2 ...
%          +kD*r_phi^3*alpha_B^(5/2)*lamt^2)*tan(phi)^2)...
%          /(4*(1+alpha_B)^2*lamt^3*besselk(1,r_phi/lambda,1));
%         parts(9) = (besselk(0,r_phi/lambda,1)^2)* ...
%          ((kD*pi*r_phi^2*tan(phi)^2*alpha_B^2)/(4*(1+alpha_B)^2) ...
%          -(kD*5*pi*r_phi^2*tan(phi)^2*alpha_B)/(8*(1+alpha_B)) ...
%          +pi*r_phi^2*kappa*tan(phi)^2/(2*lamt^2) ...
%          +(pi*r_phi^2*(4*kappa+4*kD*alpha_B*kappa-3*kD*alpha_B*lamt^2)*tan(phi)^2 ...
%          /(8*(1+alpha_B)*lamt^2))...
%          )/(besselk(1,r_phi/lambda,1))^2;
%         parts(10) = (besselk(0,r_phi/lambda,1)^3* ...
%          (pi*r_phi^3*sqrt(alpha_B)*kappa*tan(phi)^2/(4*lamt^3)  ...
%          -kD*pi*r_phi^3*alpha_B^(3/2)*tan(phi)^2/(8*(1+alpha_B)*lamt) ...
%          +pi*r_phi^3*tan(phi)^2*(2*sqrt(alpha_B)*kappa+2*alpha_B^(3/2) ...
%          *kappa-kD*alpha_B^(3/2)*lamt^2)/(8*(1+alpha_B)*lamt^3))) ...
%          /(besselk(1,r_phi/lambda,1)^3);

        alpha_B_gradient = (...
        -kD*1/2*d^2*alpha_B^2/(1+alpha_B)^2 ...
        +kD*d^2*alpha_B/(1+alpha_B) ...
        +kD/2*pi*r_phi^2*alpha_B^2/(1+alpha_B)^2 ... 
        -kD*pi*r_phi^2*alpha_B/(1+alpha_B) ...
        -kD*pi*r_phi^2*alpha_B^2*tan(phi)^2/(1+alpha_B)^2/4 ...
        +kD*3/4*r_phi^2*alpha_B*tan(phi)^2/(1+alpha_B)^2 ...
        -pi*r_phi^2*kappa*tan(phi)^2/(2*lamt^2)...
        +(pi*besselk(0,r_phi/lambda,1)* ...
        (-4*kD*r_phi*sqrt(alpha_B)*lamt^4 ...
         -2*kD*r_phi*alpha_B^(3/2)*lamt^4 ...
         -2*r_phi^3*sqrt(alpha_B)*kappa ...
         -4*r_phi^3*alpha_B^(3/2)*kappa ...
         -2*r_phi^3*alpha_B^(5/2)*kappa ...
         +kD*r_phi^3 *alpha_B^(3/2)*lamt^2 ...
         +kD*r_phi^3*alpha_B^(5/2)*lamt^2)*tan(phi)^2)...
         /(4*(1+alpha_B)^2*lamt^3*besselk(1,r_phi/lambda,1)) ...
        + (besselk(0,r_phi/lambda,1)^2)* ...
         ((kD*pi*r_phi^2*tan(phi)^2*alpha_B^2)/(4*(1+alpha_B)^2) ...
         -(kD*5*pi*r_phi^2*tan(phi)^2*alpha_B)/(8*(1+alpha_B)) ...
         +pi*r_phi^2*kappa*tan(phi)^2/(2*lamt^2) ...
         +(pi*r_phi^2*(4*kappa+4*kD*alpha_B*kappa-3*kD*alpha_B*lamt^2)*tan(phi)^2 ...
         /(8*(1+alpha_B)*lamt^2))...
         )/(besselk(1,r_phi/lambda,1))^2 ...
        + (besselk(0,r_phi/lambda,1)^3* ...
         (pi*r_phi^3*sqrt(alpha_B)*kappa*tan(phi)^2/(4*lamt^3)  ...
         -kD*pi*r_phi^3*alpha_B^(3/2)*tan(phi)^2/(8*(1+alpha_B)*lamt) ...
         +pi*r_phi^3*tan(phi)^2*(2*sqrt(alpha_B)*kappa+2*alpha_B^(3/2) ...
         *kappa-kD*alpha_B^(3/2)*lamt^2)/(8*(1+alpha_B)*lamt^3))) ...
         /(besselk(1,r_phi/lambda,1)^3));

        phi_gradient = 0;

        g = [alpha_A_gradient, alpha_B_gradient, phi_gradient]';

        E_adhesion = epsilon*n0*S_A./(1+alpha_A);
        E_stretch_A = kD/2*(alpha_A.^2*S_A./(1+alpha_A));
        E_stretch_B = kD/2*(alpha_B.^2*S_B./(1+alpha_B));
        E_bend_B = Ebend;
        E_bend_A = 4*pi*kappa*(1-cos(phi));

        % stretching, adhesion and bending energy
        f = E_adhesion + E_stretch_A + E_stretch_B + E_bend_A + E_bend_B;

    end

    function [c,ceq, gc, gceq]= constraint(~)
    
        c = [];
        ceq = S_A./(1+alpha_A)+S_B./(1+alpha_B)-S_i/(1+alpha_i);

        gc = [];

        alpha_A_gradient = -S_A/(1+alpha_A)^2;

        alpha_B_gradient = (-1).*d.^2.*(1+alpha_B).^(-2)+pi.*R.^2.*(1+alpha_B).^(-2).*sin(phi).^2+(-1/2) ...
          .*pi.*R.^2.*(1+alpha_B).^(-2).*sin(phi).^2.*tan(phi).^2+(1/2).*pi.*R.^2.* ...
          alpha_B.^(-1).*(1+alpha_B).^(-1).*sin(phi).^2.*tan(phi).^2+(-1/2).*pi.*R.^3.* ...
          alpha_B.^(-1/2).*(1+alpha_B).^(-1).*lamt.^(-1).*besselk(0,R.*alpha_B.^(1/2).*lamt.^( ...
          -1).*sin(phi),1).^3.*besselk(1,R.*alpha_B.^(1/2).*lamt.^(-1).*sin(phi),1).^(-3).* ...
          sin(phi).^3.*tan(phi).^2+(1/2).*pi.*alpha_B.^(-1/2).*(1+alpha_B).^(-2).*lamt.^(-1) ...
          .*besselk(0,R.*alpha_B.^(1/2).*lamt.^(-1).*sin(phi),1).*besselk(1,R.*alpha_B.^( ...
          1/2).*lamt.^(-1).*sin(phi),1).^(-1).*(2.*R.*lamt.^2.*sin(phi)+R.^3.*sin(phi) ...
          .^3+R.^3.*alpha_B.*sin(phi).^3).*tan(phi).^2+besselk(0,R.*alpha_B.^(1/2).*lamt.^( ...
          -1).*sin(phi),1).^2.*besselk(1,R.*alpha_B.^(1/2).*lamt.^(-1).*sin(phi),1).^(-2).* ...
          ((1/2).*pi.*R.^2.*(1+alpha_B).^(-2).*sin(phi).^2.*tan(phi).^2+(-1).*pi.* ...
          R.^2.*alpha_B.^(-1).*(1+alpha_B).^(-1).*sin(phi).^2.*tan(phi).^2);

        phi_gradient = 0;


        gceq = [alpha_A_gradient,  alpha_B_gradient, phi_gradient]';


    end

    function ceq = lipid_con_bend_test_A(alpha_B)
   
    
        % problems when phi = 0, so change slightly
        if phi==0
            phi = phi+eps(1);
        end
        
        Sigma = kD*alpha_B;
        r_phi = sin(phi)*R;
        
        [~, delA, ~] =...
            free_shape_linear_free_h(R, phi, kappa, Sigma);
    
        S_A = 2*pi*R^2*(1-cos(phi));
        S_B = delA+(d^2-pi*r_phi^2);
    
        ceq(1) = S_A./(1+alpha_A)+S_B./(1+alpha_B)-S_i/(1+alpha_i);
    
    end

end
