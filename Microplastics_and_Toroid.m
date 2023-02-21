% implements Baulin 2021 microplastics model
% all SI units
clear variables

% d_vals = linspace(1.8,1.7,4);
% for d = d_vals
% d = 1;
R = 0.3;
sigma = 5e-4;
d = sqrt(pi*R^2/sigma);
gamma = 0.3*R;            % um
alpha_i = 0.001;
kD = 300/10^12*1e9; % picoJ/um^2
kappa = 1e-19*1e12; % picoJd

% zeta = -1e-5;
% epsilon = zeta*kD; % picoJ/um^2
epsilon = -kappa/gamma^2;
zeta = epsilon/kD;

% sigma = R^2/d^2;
phi_vals = deg2rad(linspace(5,175,100));
% phi_vals = deg2rad(130);
plotfigs = 0;

if plotfigs
    figure();
    hold on
    axis equal
end
for ii=1:length(phi_vals)

    if ii==1
        alpha_A_init = alpha_i;
        alpha_B_init = alpha_i;
        rho_init = R/100;
    else
        alpha_A_init = alpha_A_vals(ii-1);
        alpha_B_init = alpha_B_vals(ii-1);
        rho_init = rho_vals(ii-1);
    end


    phi = phi_vals(ii)

    options = optimoptions('fmincon','MaxFunEvals', 1e5, 'MaxIter', 1e5, 'algorithm', 'sqp');
%     options = optimoptions('fmincon','MaxFunEvals', 1e5, 'MaxIter', 1e5,...
%         'algorithm', 'interior-point',...
%         "EnableFeasibilityMode",true);

    const = [epsilon, alpha_i, d, R, phi, kappa, kD];
    [out,fval,exitflag,output,lam_vals,grad,hessian] = ...
        fmincon(@(y) free_theta(y, const),[0.001, 0.001, R/20], [],[],[],[],...
        [-0.1,-0.1,0],[0.1, 0.1,5*d], ...
        @(y) area_con_theta(y,const), options);

    E_vals(ii) = fval;
    alpha_A_vals(ii) = out(1);
    alpha_B_vals(ii) = out(2);
    rho_vals(ii) = out(3);
    delta_vals(ii) = sin(phi)*(R+rho_vals(ii));

    if plotfigs
        t = linspace(-pi/2,-pi/2+phi,1000);
    %     t = linspace(-pi/2,pi/2,1000);
        x = cos(t)*R;
        % y = sin(t)*R+(R*cos(phi)+h(1));
        y = sin(t)*R+R*cos(phi);
        plot(x,y, 'k:', 'HandleVisibility','off')
    
        tt = linspace(-3*pi/2,-3*pi/2+phi,1000);
    %     t = linspace(-pi/2,pi/2,1000);
        xt = cos(tt)*rho_vals(ii);
        % y = sin(t)*R+(R*cos(phi)+h(1));
        yt = sin(tt)*rho_vals(ii)+rho_vals(ii)*cos(phi);
        plot(xt-xt(end)+x(end),yt-yt(end), 'k-', 'HandleVisibility','off', 'linewidth', 0.5)
        
        xf = linspace(xt(1)-xt(end)+x(end), d/2,100);
        yf = ones(size(xf))*(yt(1)-yt(end));
        plot(xf, yf, 'k-', 'HandleVisibility','off', 'linewidth', 0.5)
    end

end

%%
% anno_string = [sprintf('$k_D = %0.2g$ pJ/$\\mu$m$^2$ \n', kD),...
%     sprintf('$\\kappa = %0.2g$ pJ \n', kappa),...
%     sprintf('$\\alpha_i = %0.2g$ \n', alpha_i),...
%     sprintf('$R = %0.2g$ $\\mu$m \n', R),...
%     sprintf('$\\zeta = %0.2g$ \n', zeta)];
% 
% anno_string = [sprintf('$k_D = %0.2g$ pJ/$\\mu$m$^2$ \n', kD),...
%     sprintf('$\\kappa = %0.2g$ pJ \n', kappa),...
%     sprintf('$\\alpha_i = %0.2g$ \n', alpha_i),...
%     sprintf('$R = %0.2g$ $\\mu$m \n', R),...
%     sprintf('$\\zeta = %0.2g$ \n', zeta)];

figures = findall(groot, 'type', 'figure');

% if isempty(figures)
%     figure();
% else
%     [f1,f2] = figures.Number;
%     if f1==1
%         figure(figures(1));
%     else
%         figure(figures(2));
%     end
% end
% hold on
% plot(phi_vals, (E_vals-E_vals(1))/(kD));
% % ylim([min(E_vals), max(E_vals)])
% xlabel('$\phi$')
% ylabel('$E_\mathrm{total}$')

% relative shapes
if isempty(figures)
    figure();
else
    [f1,f2] = figures.Number;
    if f1==1
        figure(figures(1));
    else
        figure(figures(2));
    end
end
hold on
% plot(phi_vals, (E_vals-E_vals(1))/((E_vals(end)-E_vals(1))),...
%     'displayname', sprintf('$\\sigma = %0.2g$', sigma));
plot(phi_vals, (E_vals-E_vals(1))/max(abs(E_vals)),...
    'displayname', sprintf('$\\alpha_i = %0.2g$', alpha_i));
% ylim([min(E_vals), max(E_vals)])
xlabel('$\phi$')
ylabel('$E_\mathrm{total}$')
legend

newcolors = ["#016D24","#648FFF","#FE6100","#DC267F","#1D0C6B","#27E0D8","#FFB000"];
colororder(newcolors)

% legend(sprintf('$\\alpha_i = %d, \\alpha_B = %d$', ))

% annotation('textbox', [0.3,0.6,0.4,0.2], 'String',...
%     anno_string,...
%     'interpreter', 'latex','FontSize',16,...
%     'FitBoxToText','on','LineStyle','none', 'Color','b')

% figure();
% hold on
% plot(rad2deg(phi_vals), rho_vals)

if isempty(figures)
    figure();
else
    if f1==2
        figure(figures(1));
    else
        figure(figures(2));
    end
end
hold on
plot(rad2deg(phi_vals), alpha_B_vals)
xlabel('$\phi$')
ylabel('$\alpha_B$')

% end

newcolors = ["#016D24","#648FFF","#FE6100","#DC267F","#1D0C6B","#27E0D8","#FFB000"];
colororder(newcolors)


%% functions

function f = free_theta(y, const)
    % function of free energy to minimise

    epsilon = const(1);
    alpha_i = const(2);
    d = const(3);
    R = const(4);
    phi = const(5);
    kappa = const(6);
    kD = const(7);

    alpha_A = y(1);
    alpha_B = y(2);
    rho = y(3); %toroid radius

    delta = sin(phi)*(R+rho);

    S_A = 2*pi*R^2*(1-cos(phi));
    S_B = d^2 - pi*delta^2 + 2*pi*rho*(delta*phi-rho*(1-cos(phi)));

    E_adhesion = epsilon*S_A./(1+alpha_A);
    E_stretch_A = kD/2*alpha_A.^2*S_A./(1+alpha_A);
    E_stretch_B = kD/2*alpha_B.^2*S_B./(1+alpha_B);
    E_bend_A = 4*pi*kappa*(1-cos(phi));
    a = delta/rho;
    if a>1
        E_bend_B = 2*pi*kappa*(a^2/sqrt(a^2-1)*...
            (acot(sqrt(a^2-1))+atan((a*tan(phi/2)-1)/sqrt(a^2-1)))-2*(1-cos(phi)));
    else
        E_bend_B = 2*pi*kappa*(a^2/sqrt(1-a^2)*...
            (-acoth(sqrt(1-a^2))+atanh((1-a*tan(phi/2))/sqrt(1-a^2)))-2*(1-cos(phi)));
    end

    f = E_adhesion+E_stretch_A+E_stretch_B+E_bend_A+E_bend_B;

end

function [c,ceq] = area_con_theta(y, const)

    epsilon = const(1);
    alpha_i = const(2);
    d = const(3);
    R = const(4);
    phi = const(5);
    kappa = const(6);
    kD = const(7);

    alpha_A = y(1);
    alpha_B = y(2);
    rho = y(3);

    delta = sin(phi)*(R+rho);

    S_i = d^2;
    S_A = 2*pi*R^2*(1-cos(phi));
    S_B = d^2 - pi*delta^2 + 2*pi*rho*(delta*phi-rho*(1-cos(phi)));

    c(1) = 0;

    ceq(1) = S_A./(1+alpha_A)+S_B./(1+alpha_B)-S_i/(1+alpha_i);

end