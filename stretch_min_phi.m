function f = stretch_min_phi(y, const)
    % function of free energy to minimise

    zeta = const(1);
    alpha_i = const(2);
    d = const(3);
    R = const(4);
    phi = const(5);
    S_B = const(6);

    alpha_A = y(1);
    alpha_B = y(2);

    S_A = 2*pi*R^2*(1-cos(phi));

    f = zeta*S_A./(1+alpha_A) ...
      + 1/2*(alpha_A.^2*S_A./(1+alpha_A)+alpha_B.^2*S_B./(1+alpha_B));

end