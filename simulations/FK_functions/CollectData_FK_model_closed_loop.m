function [X,Y,U] = CollectData_FK_model_closed_loop(Ntraj, Tsim, Ts, f_u, X_ref, X0, fk_param_struct)
% --- Inputs ---
% Ntraj:    number of trajectories
% Tsim:     Simulation time of a single trajectory
% Ts:       Sampling period
% f_u:      Closed-loop controller, function handler
% X_ref:    Reference position for all pendulums for a single trajectory
% X0:       Initial positions for all pendulums
% fk_param_struct

% --- Outputs ---
% X: All trajectories points (states)
% Y: All trajectories points related to the X 
% U: Inputs related to the X

%%% Notes:
% This function is not applicable for non-constant trajectory
% Assume input of size 1

N = fk_param_struct.N;
m = fk_param_struct.m;
l = fk_param_struct.l;
gamma = fk_param_struct.gamma;
b = fk_param_struct.b;
k = fk_param_struct.k;
J = fk_param_struct.J;
L = fk_param_struct.L;
d = fk_param_struct.d;
g = fk_param_struct.g;
f_origin = fk_param_struct.f_origin;

Nsim = Tsim/Ts;
t = 0:Ts:Tsim;  % time vector for a single trajectory

% Initialize All matrices for all trajectories
X = zeros(size(X0,1),Ntraj*Nsim);
Y = zeros(size(X0,1),Ntraj*Nsim);
U = zeros(1,Ntraj*Nsim);

% For every trajectory
for traj = 1:Ntraj
    fprintf('Trajectory %d out of %d \n',traj,Ntraj);
    
    x0 = X0(:,traj);
    Xsim = zeros(size(x0,1), Nsim+1);
    Xsim(:,1) = x0;
    Usim = zeros(1, Nsim);
    % and now the closed-loop simulation
    for ii = 2:Nsim+1
        tspan = [t(ii-1), t(ii)];
        u = f_u(x0, X_ref(:,ii));
        u1 = u(1,:); % Select only the first step of the input (from the MPC) 
        u_disturbance = zeros(N,1); % No external disturbances to pendulums during identification as it SHOULD BE
        [~, Y_curr] = ode45(@(t,x) FK_solver(t, x, u1, u_disturbance, N, L, d, f_origin, m, l, gamma, b, k, J, g), tspan, x0);

        Xsim(:,ii) = Y_curr(end,:)';
        Usim(ii-1) = u(1,:);
        x0 = Y_curr(end,:)';
    end
    
%     % Debug:
%     figure;
%     plot(Xsim(1:2:end, :)');
%     hold on;
%     plot(X_ref(1,:), 'Linewidth', 2);
%     plot(Usim);
    
    % Store data
    idx_start = (traj-1)*Nsim + 1;
    idx_end = (traj)*Nsim;
    U(:,idx_start:idx_end) = Usim(1:end);
    X(:,idx_start:idx_end) = Xsim(:,1:end-1);
    Y(:,idx_start:idx_end) = Xsim(:,2:end);
end







