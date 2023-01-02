clear variables

R = 0.3;
alpha_i = 0.01;
st = alpha_i;
d = R*10;

ptmatrix = [-st, -st, 5, 0;-st, st, 5, 0;st, -st, 5, 0;st, st, 5, 0;...
    -st, -st, 5, R;-st, -st, 5, -R;st, st, 5, d/2;st, st, 5, -d/2;...
    -st, -st, 85, 0;-st, st, 85, 0;st, -st, 85, 0;st, st, 85, 0;...
    -st, -st, 85, R;-st, -st, 85, -R;st, st, 85, d/2;st, st, 85, -d/2;...
    -st, -st, 45, 0;-st, st, 45, 0;st, -st, 45, 0;st, st, 45, 0;...
    -st, -st, 45, R;-st, -st, 45, -R;st, st, 45, d/2;st, st, 45, -d/2];


alpha_A_vals_inp = [-0.1, -alpha_i,0, alpha_i, 0.1];
alpha_B_vals_inp = [-0.1, -alpha_i,0, alpha_i, 0.1];
phi_vals_inp = [0.001, 5, 45, 85, pi/2-0.001];
h_phi_vals_inp = [d/2, R, 0, -R, -d/2];
for aa = 1:length(alpha_A_vals_inp)
    for bb = 1:length(alpha_B_vals_inp) 
        for pp = 1:length(phi_vals_inp)
            for hh = 1:length(h_phi_vals_inp) 
                ptmatrix_test(aa,bb,pp,hh,1:4) = ...
                    [alpha_A_vals_inp(aa), alpha_A_vals_inp(bb),...
                     phi_vals_inp(pp), h_phi_vals_inp(hh)];
            end
        end
    end
end

