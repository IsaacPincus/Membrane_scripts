function fig_handle = dimensionless_parameters(xaxis_name, xscale, yaxis_name, yscale,...
    values_structure)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

param_1_vals = values_structure.param_1_vals';
lambda_vals = values_structure.lambda_vals;
sigma_vals = values_structure.sigma_vals;
d_vals = values_structure.d_vals;

fig_handle = figure('Position',[400,100,800,600]);
hold on
axes1 = gca;
axes1.XScale = xscale;
axes1.YScale = yscale;
plot(param_1_vals, param_1_vals./lambda_vals, 'DisplayName','$R/\lambda$')
plot(param_1_vals, sigma_vals, 'DisplayName','$\sigma$')
plot(param_1_vals, d_vals./lambda_vals, 'DisplayName','$d/\lambda$')
xlabel(xaxis_name)
ylabel(yaxis_name)
legend
end