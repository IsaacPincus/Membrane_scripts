N = 1e3;
d = 2;
R = 0.2;
phi = deg2rad(30);
r_phi = sin(phi)*R;
h_phi = -0.1;
r = linspace(r_phi, d/2,N);

figure();
hold on
xlabel('$r$')
ylabel('$h$')
% ylabel('$\partial h/\partial r$')
% ylabel('$\nabla^2 h$')
% axis equal
colours = ["r", "g", "b", "m"];
ii = 0;

for lambda=[1,0.1,0.01,0.005]
% for lambda=[0.01]
    ii = ii+1;
    
    % full solution
    A_c = d/(2*r_phi*lambda)*besselk(1,d/2/lambda) - 1/lambda*besselk(1,r_phi/lambda);
    B_c = -d/(2*r_phi*lambda)*besseli(1,d/2/lambda) + 1/lambda*besseli(1,r_phi/lambda);
    D_c = besseli(0,r_phi/lambda) - besseli(0,d/2/lambda) ...
        - log(2*r_phi/d)*d/2/lambda*besseli(1,d/2/lambda);
    E_c = besselk(0,r_phi/lambda) - besselk(0,d/2/lambda) ...
        + log(2*r_phi/d)*d/2/lambda*besselk(1,d/2/lambda);
    
    C(4) = (h_phi - D_c*tan(phi)/B_c)/(-A_c*D_c/B_c+E_c);
    C(3) = (tan(phi) - C(4)*A_c)/B_c;
    C(1) = -C(3)*d/2/lambda*besseli(1,d/2/lambda) + C(4)*d/2/lambda*besselk(1,d/2/lambda);
    C(2) = -C(1)*log(d/2/lambda) - C(3)*besseli(0,d/2/lambda) - C(4)*besselk(0, d/2/lambda);
    
    h = C(1)*log(r/lambda) + C(2) + C(3)*besseli(0,r/lambda) + C(4)*besselk(0,r/lambda);

    hderiv = C(1)./r+C(3)/lambda*besseli(1,r/lambda)-C(4)/lambda*besselk(1,r/lambda);
    lap_h = C(3)/lambda^2*besseli(0,r/lambda)+C(4)/lambda^2*besselk(0,r/lambda);
    
    plot(r,h, strcat(colours(ii),'-'), ...
        'displayname', sprintf('$\\lambda = %0.2g$', lambda));

%     plot(r,hderiv, strcat(colours(ii),'-'), ...
%         'displayname', sprintf('$\\lambda = %0.2g$', lambda));

%     plot(r,lap_h, strcat(colours(ii),'-'), ...
%         'displayname', sprintf('$\\lambda = %0.2g$', lambda));

    % asymptotic solution
    C4_bar = (d/(2*r_phi*lambda)*besseli(1,d/2/lambda,1)*h_phi ...
        - besseli(0,d/2/lambda,1)*tan(phi) -...
        log(2*r_phi/d)*d/(2*lambda)*besseli(1,d/2/lambda,1)*tan(phi))/...
        (d/(2*r_phi*lambda)*besselk(0,r_phi/lambda,1)*besseli(1,d/2/lambda,1)...
        + 1/lambda*besselk(1,r_phi/lambda,1)*besseli(0,d/2/lambda,1)...
        + log(2*r_phi/d)*d/(2*lambda^2)*besselk(1,r_phi/lambda,1)*besseli(1,d/2/lambda,1));
    C3_bar = (tan(phi)+1/lambda*besselk(1,r_phi/lambda,1)*C4_bar)/...
        (-d/(2*r_phi*lambda)*besseli(1,d/2/lambda,1));
    h_as = -d/2/lambda*C3_bar*besseli(1,d/2/lambda,1)*log(r/lambda)...
        +d/2/lambda*C3_bar*log(d/2/lambda)*besseli(1,d/2/lambda,1)...
        -C3_bar*besseli(0,d/2/lambda,1)...
        +exp((r-d/2)/lambda).*C3_bar.*besseli(0,r/lambda,1)...
        +exp((r_phi-r)/lambda).*C4_bar.*besselk(0,r/lambda,1);

    hderiv_as = -besseli(1,d/2/lambda,1)*C3_bar*d/(2*lambda)*1./r...
        +exp((r-d/2)/lambda).*C3_bar/lambda.*besseli(1,r/lambda,1)...
        -exp((r_phi-r)/lambda).*C4_bar/lambda.*besselk(1,r/lambda,1);
    lap_h_as = exp((r-d/2)/lambda).*C3_bar/lambda^2.*besseli(0,r/lambda,1)...
        +exp((r_phi-r)/lambda).*C4_bar/lambda^2.*besselk(0,r/lambda,1);

    plot(r,h_as, strcat(colours(ii),':'), ...
        'displayname', sprintf('$\\lambda = %0.2g$ asymptotic', lambda));

%     plot(r,hderiv_as, strcat(colours(ii),':'), ...
%         'displayname', sprintf('$\\lambda = %0.2g$ asymptotic', lambda));

%     plot(r,lap_h_as, strcat(colours(ii),':'), ...
%         'displayname', sprintf('$\\lambda = %0.2g$ asymptotic', lambda));
end

% laplace solution
hl = h_phi*log(2*r/d)/log(2*r_phi/d);

% plot(r,hl,'k-.', 'displayname','$\lambda \rightarrow 0$');

legend






















