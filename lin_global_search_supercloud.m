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
R_vals = 0.3;                                   % um
sigma = 0.01;                                   % surface fraction
d_vals = sqrt(R_vals.^2/sigma);                 % um
kD_vals = 300/10^12*1e9;                        % picoJ/um^2
zeta = 0.02;                                    % dimensionless
epsilon_vals = -zeta*kD_vals;                        % picoJ/um^2
n0_vals = 1;                                    % fraction
kappa_vals = logspace(-21,-15, 12)*1e12;         % picoJ
% kappa_vals = logspace(-21,-15, 60)*1e12;
alpha_i_vals = 0.01;                            % fraction
% other constants
N = 3e3;

ii = 0;
for rr = 1:length(R_vals)
    for dd = 1:length(d_vals)
        for kk = 1:length(kD_vals)
            for ee = 1:length(epsilon_vals)
                for nn = 1:length(n0_vals)
                    for pp = 1:length(kappa_vals)
                        for aa = 1:length(alpha_i_vals)
                            ii = ii+1;
                            parameter_set(ii, 1:7) = [epsilon_vals(ee),...
                                n0_vals(nn),...
                                d_vals(dd),...
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


% split up job between task IDs
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
    h_phi_init = 0;
    S_A_init = 2*pi*R^2*(1-cos(phi_init));
    S_B_init = d^2 - pi*R^2*sin(phi_init)^2;

    E_adhesion_init = epsilon*n0*S_A_init ./(1+alpha_A_init );
    E_stretch_A_init  = kD/2*(alpha_A_init .^2*S_A_init ./(1+alpha_A_init ));
    E_stretch_B_init  = kD/2*(alpha_B_init .^2*S_B_init ./(1+alpha_B_init ));
    E_total_init = E_adhesion_init+E_stretch_A_init+E_stretch_B_init;
    
    % use initial stretch to solve for minimum attachment
    const = [epsilon, n0, d, R, kD, kappa, alpha_i, N];
    opts = optimoptions(@fmincon,'Algorithm',"sqp");
%     rng default % For reproducibility
    problem = createOptimProblem("fmincon",...
        'x0', [alpha_A_init, alpha_B_init, phi_init, h_phi_init],...
        'objective', @(y) stretch_bend_min(y, const),...
        'lb', [-0.1,-0.1,0,-2*d],...
        'ub', [0.1, 0.1, pi/2, 2*d],...
        'nonlcon', @(y) lipid_con_bend(y,const),...
        'options', opts);
    alpha_A_vals_inp = [-0.1, -alpha_i, -1e-3, 0, 1e-3, alpha_i, 0.1];
%     alpha_B_vals_inp = [-0.1, -alpha_i, -1e-3, 0, 1e-3, alpha_i, 0.1];
    phi_vals_inp = deg2rad([0, 5, 20, 45, 60, 85, 90]);
    h_phi_vals_inp = [-2*d,-d/2, -R, 0, R, d/2,2*d];
    ptmatrix_this_run =...
            zeros(length(alpha_A_vals_inp),length(phi_vals_inp),...
                  length(h_phi_vals_inp), 4);
    tic
    for aa = 1:length(alpha_A_vals_inp)
        for pp = 1:length(phi_vals_inp)
            for hh = 1:length(h_phi_vals_inp) 
                const_con = [epsilon, n0, d, R, kD, kappa, alpha_i, N,...
                    alpha_A_vals_inp(aa), phi_vals_inp(pp), h_phi_vals_inp(hh)];
                alpha_B_vals_inp(pp,hh,aa) =...
                    fzero(@(x) lipid_con_bend_test_A(x, const_con), rand()*0.2-0.1);
                ptmatrix_this_run(aa,pp,hh,1:4) = ...
                    [alpha_A_vals_inp(aa), alpha_B_vals_inp(pp,hh,aa),...
                     phi_vals_inp(pp), h_phi_vals_inp(hh)];
            end
        end
    end
    toc
    ptmatrix = reshape(ptmatrix_this_run, [numel(ptmatrix_this_run)/4, 4]);
    tpoints = CustomStartPointSet(ptmatrix);
    rs = RandomStartPointSet('NumStartPoints',100);
    gs = MultiStart("FunctionTolerance",1e-3, "XTolerance", 1e-3);
    tic
    [out,fval2,exitflag,output,solutions] = run(gs,problem, {tpoints,rs});
    toc

    solution_set{ii} = solutions;
    
    % get energies etc
    alpha_A = out(1);
    alpha_B = out(2);
    phi = out(3);
    h_phi = out(4);
    
    % get the shape of the free region
    Sigma = kD*alpha_B;
    lambda = sqrt(kappa/Sigma);
    r_phi = sin(phi)*R;
    r = linspace(r_phi, d/2,N);
    
    [h,C,S_B, ~, lap_h, hderiv] = free_shape_linear_fixed_h(r, r_phi, d, phi, kappa, Sigma, h_phi);
    
    S_A = 2*pi*R^2*(1-cos(phi));
    
    E_adhesion = epsilon*n0*S_A./(1+alpha_A);
    E_stretch_A = kD/2*(alpha_A.^2*S_A./(1+alpha_A));
    E_stretch_B = kD/2*(alpha_B.^2*S_B./(1+alpha_B));
    E_bend_B = kappa/2*2*pi*trapz(r, r.*lap_h.^2);
    E_bend_A = 4*pi*kappa*(1-cos(phi));
    
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
    h_phi_vals(ii) = h_phi;
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

function ceq = lipid_con_bend_test_A(alpha_B, const)

    epsilon = const(1);
    n0 = const(2);
    d = const(3);
    R = const(4);
    kD = const(5);
    kappa = const(6);
    alpha_i = const(7);
    N = const(8);
    alpha_A = const(9);
    phi = const(10);
    h_phi = const(11);

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

    ceq(1) = S_A./(1+alpha_A)+S_B./(1+alpha_B)-S_i/(1+alpha_i);

end
