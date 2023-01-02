function dydt = FK_solver(t, y, u, L, D, f_origin, m, l, gamma, b, k, J, g)
% u: Torque acting on the first pendulum in the chain, numerical value
% Origin : 0: downward, 1: upward

x = y(1:2:end);
v = y(2:2:end);

A_inter = kron(-L, [0, 0; k, b])/J;

if f_origin == 0
    f_drift = kron(v,[1; 0]) - kron(g/l*sin(x) + (gamma/J)*v, [0; 1]);
else 
    f_drift = kron(v,[1; 0]) - kron( -g/l*sin(x) + (gamma/J)*v, [0; 1]);
end

dydt = f_drift + A_inter*y + kron(diag(D), [0;1])/J*u; % 

end

