function [T,Y] = build_and_sim_nonlin_FK_model(fk_params, x_curr, u, u_disturbance, tspan)

% With disturbance

N = fk_params.N;
m = fk_params.m;
l = fk_params.l;
gamma = fk_params.gamma;
b = fk_params.b;
k = fk_params.k;
J = fk_params.J;
L = fk_params.L;
d = fk_params.d;
g = fk_params.g;
f_origin = fk_params.f_origin;

[T, Y] = ode45(@(t,x) FK_solver(t, x, u, u_disturbance, N, L, d, f_origin, m, l, gamma, b, k, J, g), tspan, x_curr);


end

