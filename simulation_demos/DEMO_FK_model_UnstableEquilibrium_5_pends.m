%% DEMO: synchronization to unstable equilibrium: swing up
% author: Loi Do
% doloi@fel.cvut.cz

% Demo, showing synchronization to unstable equilibrium of 5 pendulums

clc;
clear;
rng(10); % Select a seed for random numbers

%% ************************ Parameters ********************** %%
%***** USER-DEFINED TUNING PARAMETERS *****

%%% FK model
N = 5; % Number of pendulums
fk_params.f_origin = 0; % Boolean flag: origin of the system in 0 (downward position); 1 (upward position); 
fk_params.N = N;  % Number of pendulums
fk_params.m = 0.017;                            % Array of weigths of a pendulum [kg];
fk_params.l = 0.15;                             % Length of a pendulum
fk_params.gamma = 3.195433320423741e-04;        % Absolute friction coeficient 
fk_params.b = 6.819092300095237e-04;            % Relative friction coeficient 
fk_params.k = 0.047966777081821;                % Spring stiffness coeficient 
fk_params.J = fk_params.m*fk_params.l^2;        % Moment of inertia of a pendulum
fk_params.g = 9.81;                             % Gravity constant

%%% MISC
Ts = 0.005;     % Sampling time for discrete-time system [seconds]

%%% Identification
%%%% LQR
Q_LQR = kron(eye(N), [1000 0; 0 0.01]); % Penalization of system's states
R_LQR = 0.1;
u_LQR_random_sigma = 0.1; % Amplitude of random perturbation to 

Ntraj = 200;    % Number of trajectories for identification
Tsim_idf = 1;    % Simulation time of a single trajectory for identification
X0_idf =  repmat([0;0],N,1) + 2*pi*repmat([1;0],N,1).*(rand(2*N,Ntraj));    % Initial states for all trajectories with random perturbations
Xref = repmat(kron(ones(N,1), [pi;0]), 1, Tsim_idf/Ts + 1); % Reference for all pendulums and for whole time interval

%%% MPC
Qy = kron(eye(N), diag([0,0.01])); % Penalization of system's outputs
idx = 1;
for ii=1:2:2*N
    Qy(ii,ii) = 10*idx^3;
    idx = idx + 1;
end
R = 0.1;

umin =-0.1; % Lower bounds on control
umax = 0.1; % Upper bounds on control
fk_params.umin = umin;
fk_params.umax = umax;

Npred = 50; % Prediction horizon´for MPC
Tsim = 3; % Simulation length of MPC 

%***** COMPUTED PARAMETERS *****
adj_mat = (diag(ones(1,N-1), 1) + diag(ones(1,N-1), -1)); % Adjacency matrix
L = (diag(sum(adj_mat)) - adj_mat);
D = zeros(size(L));
D(1,1) = 1;
fk_params.L = L; % System's laplacian
fk_params.D = D; % System's input description

%% ************************PART I: Identification ********************** %%
% Create controller based on linearized FK model.
%%% Linearized system around system's origin 
[~, sys_c] = create_linearized_FK_model(fk_params, Ts);

%%% Build LQR controller
[K_lqr,~,~] = lqr(sys_c.A,sys_c.B,Q_LQR,R_LQR,0); % Computes u=-K_lqr*x
f_u_LQR = @(x) K_lqr*x + normrnd(0, u_LQR_random_sigma); % LQR state feedback with random perturbation.
% TODO: Validate closed-loop control on a single trajectory

%%% Collecting closed-loop data for EDMD
[X,Y,U] = CollectData_FK_model_closed_loop(Ntraj, Tsim_idf, Ts, f_u_LQR, Xref, X0_idf, fk_params);
f_name = ['./data/DEMO_data_swingUp_numTraj-' , num2str(Ntraj), '.mat'];
disp('Data collected');

% EDMD
[A,B,C,f_BuildKoopmanState]= SystemID_via_EDMD_FK(X,Y,U);
%% MPC ctrl
% MPC setup for the Koopman linear system
% z'=Az+Bu ; y = Cz

% ---------------- SETUP ----------------
n = size(A,1);  % state dimension
r = size(C,1);  % output dimension
nu = size(B,2); % input dimension

% Input constraints
umin_mpc = umin.*ones(nu,1);
umax_mpc = umax.*ones(nu,1);

xlift_min = -1000*ones(n,1); % Lower bounds on states
xlift_max = 1000*ones(n,1); % Upper bounds on states

[~,~,kmpc]  = osqp_MPC_controller(A,B,C,Qy,R,Qy,Npred,umin_mpc, umax_mpc, xlift_min, xlift_max);

Nsim = Tsim/Ts;
t = 0:Ts:Tsim;

Xsim = zeros(2*N, Nsim+1); % States System starting at the origin
Usim = zeros(1, Nsim); % Inputs

% and now the closed-loop simulation
OutputError=[];
for i = 1:Nsim
    if(mod(i,10) == 0)
        fprintf('Closed-loop simulation, %f %% completed \n', 100 * i / Nsim);
    end
    
    yref = Xref(:,1);    % reference output without preview
    z = f_BuildKoopmanState(Xsim(:,i));   % create the state of the Koopman linear system via embedding and lift
    
    [u, ~] = kmpc(z,yref); % get the control input
    [~,Y_nonlin] = build_and_sim_nonlin_FK_model(fk_params, Xsim(:,i), u(1), [0 Ts]); % advance in time for one step 
    
    % store data
    Xsim(:,i+1) = Y_nonlin(end,:)';
    Usim(i) = u(1);
    
end

%% Display results
figure;
title('Synchronization to unstable equilibrium: Swing up');
subplot(2,1,1);
hold on;
idx = 1;
for ii = 1:2:2*N
    plot(t, Xsim(ii,:)', 'Linewidth', 1.5, 'Color', color_p(idx));
    idx = idx + 1;
end
plot(t, ones(size(t))*pi,'Linewidth', 3, 'Linestyle', ':');
legend('#1', '#2', '#3','#4','#5', 'Reference');
box on;
grid on;
ylabel('Angle [rad]');
xlabel('Time [s]');

subplot(2,1,2);
plot(t(1:end-1), Usim, 'Linewidth', 1.5);
legend('Input');
box on;
grid on;
ylabel('Torque [N.m]');
xlabel('Time [s]');


