function [T,Y] = build_and_sim_nonlin_FK_model(fk_params, x_curr, input, tspan)

m = fk_params.m;
l = fk_params.l;
gamma = fk_params.gamma;
b = fk_params.b;
k = fk_params.k;
J = fk_params.J;
L = fk_params.L;
D = fk_params.D;
g = fk_params.g;
f_origin = fk_params.f_origin;

options = odeset('RelTol',1e-10, 'AbsTol', 1e-10);

[T, Y] = ode45(@(t,x) FK_solver(t, x, input, L, D, f_origin, m, l, gamma, b, k, J, g), tspan, x_curr, options);


end

