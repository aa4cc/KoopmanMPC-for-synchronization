function dxdt = FK_solver(t, y, u, u_disturbance, N, L, d, f_origin, m, l, gamma, b, k, J, g)
% t: current time
% y: current state
% u: Torque acting on the first pendulum in the chain
% u_disturbance: Array of torques acting on each pendulum in the system
% N: number of pendulums
% L: graph laplacian defining the system
% d: array defining the input dynamics
% f_origin: flag: origin of the system : 0: downward, 1: upward
% m: mass of a single pendulum
% l: length of a single pendulum
% gamma: damping coefficient of the absolute speed
% b: dampáing coefficient of the relative speed
% k: spring stiffness
% J: moment of inertia of a pendulum
% g: gravity constant

x = y(1:2:end); % angles
v = y(2:2:end); % angular speeds

% Input matrices in the Multi-agent system description
G = [0;1];
K = [k, b];

% Drift (uncoupled) dynamics
if f_origin == 0
    f_drift = kron(v,[1; 0]) + kron(-m*g*l/J*sin(x) - (gamma/J)*v, [0; 1]);
else 
    f_drift = kron(v,[1; 0]) + kron(+m*g*l/J*sin(x) - (gamma/J)*v, [0; 1]);
end

A_inter = -kron(L, G*K)/J; % Linear interconnection between the pendulums
controlled_input = kron(d, G)/J*u;
disturbance_input = kron(eye(N), G)/J*u_disturbance;

dxdt = f_drift + A_inter*y + controlled_input + disturbance_input;

end

