function [X,Y,U] = CollectData_FK_model_closed_loop(Ntraj, Tsim, Ts, f_u, X_ref, X0, fk_param_struct)
% Ntraj:        number of trajectories
% N:            number of pendulums
% Tsim:         Simulation time
% Ts:           Sampling period
% f_u:          Closed-loop controller, function handler
% x_ref:        Reference position for a single pendulum
% x_0_avg:      Initial position for a single pendulum

% X: States in the current time
% Y: States in the next time
% U: Inputs in the current time

% Initialize All matrices for all trajectories
X = [];
Y = [];
U =[];

N = fk_param_struct.N;
m = fk_param_struct.m;
l = fk_param_struct.l;
gamma = fk_param_struct.gamma;
b = fk_param_struct.b;
k = fk_param_struct.k;
J = fk_param_struct.J;
L = fk_param_struct.L;
D = fk_param_struct.D;
g = fk_param_struct.g;
f_origin = fk_param_struct.f_origin;
umin = fk_param_struct.umin;
umax = fk_param_struct.umax;

Nsim = Tsim/Ts;
t = 0:Ts:Tsim;  % time vector for a single trajectory

% TODO: make this applicable for non-constant trajectory
% X_ref_stacked = repmat(x_ref, N, 1); % Reference trajectory for all pendulums

% For every trajectory
for traj = 1:Ntraj
    fprintf('Trajectory %d out of %d \n',traj,Ntraj);
    
    x0 = X0(:,traj);
    Xsim = x0;
    Usim = [];
    % and now the closed-loop simulation
    for ii = 2:Nsim+1
        tspan = [t(ii-1), t(ii)];
        u = min(max(f_u(X_ref(:,ii)-x0), umin), umax);    % Saturated input;
        [~, Y_curr] = ode45(@(t,x) FK_solver(t, x, u, L, D, f_origin, m, l, gamma, b, k, J, g), tspan, x0);

        Xsim = [Xsim, Y_curr(end,:)'];
        Usim = [Usim, u];
        x0 = Y_curr(end,:)';
    end
    
%     % Debug:
%     figure;
%     plot(Xsim(1:2:end, :)');
%     hold on;
%     plot(X_ref(1,:), 'Linewidth', 2);
%     plot(Usim)
%     
    % Store
    U = [U Usim(1:end)]; % Take from time 0 to T-dt
    X = [X Xsim(:,1:end-1)]; % Take from time 0 to T-dt
    Y = [Y Xsim(:,2:end)]; % Take from time 0+dt to T
    
end







