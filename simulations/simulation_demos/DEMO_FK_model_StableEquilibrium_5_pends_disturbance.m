%% DEMO: synchronization to stable equilibrium with disturbances
% author: Loi Do
% doloi@fel.cvut.cz

% Demo, showing synchronization to stable equilibrium of FK model with 5 
% pendulums and with external disturbances

clc;
clear;
rng(10); % Select a seed for random numbers

%% ************************ Parameters ********************** %%
%%% FK model
N = 5; % Number of pendulums
fk_params.f_origin = 0; % Boolean flag: origin of the system in 0 (downward position); 1 (upward position);                                
fk_params.m = 0.017;                            % Array of weigths of a pendulum [kg];
fk_params.l = 0.15;                             % Length of a pendulum
fk_params.gamma = 3.195433320423741e-04;        % Absolute friction coeficient 
fk_params.b = 6.819092300095237e-04;            % Relative friction coeficient 
fk_params.k = 0.047966777081821;                % Spring stiffness coeficient 
fk_params.J = fk_params.m*fk_params.l^2;        % Moment of inertia of a pendulum
fk_params.g = 9.81;                             % Gravity constant

adj_mat = (diag(ones(1,N-1), 1) + diag(ones(1,N-1), -1)); % Adjacency matrix
L = (diag(sum(adj_mat)) - adj_mat);
fk_params.L = L; % System's laplacian
fk_params.d = [1; zeros(N-1,1)]; % System's controller input description in Multi-agent system formalism
fk_params.N = N; 

%%% MISC
Ts = 0.005;     % Sampling time for discrete-time system [seconds]

%%% MPC / KMPC
umin =-0.1; % Lower bounds on control
umax = 0.1; % Upper bounds on control
fk_params.umin = umin;
fk_params.umax = umax;

Qy = kron(eye(N), diag([0,0.01])); % Penalization of system's outputs in MPC for identification
idx = 1;
for ii=1:2:2*N
    Qy(ii,ii) = 10*idx^3;
    idx = idx + 1;
end
R = 0.1;

Npred = 3; % Prediction horizon for MPC and KMPC

nu = size(R,2); % input dimension
umin_mpc = umin.*ones(nu,1);
umax_mpc = umax.*ones(nu,1);

%%% Identification
Ntraj_idf = 100;    % Number of trajectories for identification
Tsim_idf = 1;    % Simulation time of a single trajectory for identification
X0_idf =  repmat([0;0],N,1) + pi*repmat([1;0],N,1).*(rand(2*N,Ntraj_idf) - 0.5);    % Initial states for all trajectories with random perturbations
Xref_idf = repmat(kron(ones(N,1), [0;0]), 1, Tsim_idf/Ts + 1); % Reference for all pendulums and for whole time interval of identification
u_MPC_random_sigma = 0.01;  % Perturbation to the input used for identification

x_min = -1000*ones(2*N,1); % Lower bounds on states
x_max = 1000*ones(2*N,1); % Upper bounds on states

%%% Final simulation
n = 6*N;  % lifted state dimension
xlift_min = -1000*ones(n,1); % Lower bounds on lifted states
xlift_max = 1000*ones(n,1); % Upper bounds on lifted states

Tsim = 10; % Simulation length of the final simulation
Nsim = Tsim/Ts;
t = 0:Ts:Tsim;

% Disturbance for KMPC experiment: 
% ------ NOTE: change, if different number ofpendulums is used ------
disturbance_idx = [0; 2; 4; 6; 8]/Ts + 1;   % Simulation time-index of disturbance occurance
disturbance_amp = [-2; 1; 1.5; -0.5; 2];    % Amplitude of the disturbance (N.m)
disturbance_pend_idx = [1; 3; 5; 2; 4];     % Index of the pendulum being disturbed


%% ************************PART I: Identification ********************** %%
% Create controller based on linearized FK model.
[sys_d, ~] = create_linearized_FK_model(fk_params, Ts); % Linearized system around system's origin 

% Build MPC controller with linearized system
[~,~,kmpc]  = osqp_MPC_controller(sys_d.A,sys_d.B,sys_d.C,Qy,R,Qy,Npred,umin_mpc, umax_mpc, x_min, x_max);
f_kmpc = @(x, r) kmpc(x,r) + normrnd(0, u_MPC_random_sigma);    % Add random perturbation

% Collecting closed-loop data for EDMD
f_collect_data = 0; % Flag: 1: collect data, 0: load saved data
if f_collect_data == 1
    [X,Y,U] = CollectData_FK_model_closed_loop(Ntraj_idf, Tsim_idf, Ts, f_kmpc, Xref_idf, X0_idf, fk_params);
    f_name = ['./simulations/data/DEMO_data_StableEq_5_pends_numTraj-' , num2str(Ntraj_idf), '.mat'];
    save(f_name, 'X', 'Y', 'U');
else
    load('./simulations/data/DEMO_data_stableEq_MPC_idf_numTraj-100.mat');  
end
disp('Data collected');

% Run EDMD
[A,B,C,f_BuildKoopmanState]= SystemID_via_EDMD_FK(X,Y,U); 
%% MPC ctrl
% MPC setup for the Koopman linear system
% z_{k+1}=Az_{k}+Bu_{k} ; y_{k} = Cz_{k}

[~,~,kmpc]  = osqp_MPC_controller(A,B,C,Qy,R,Qy,Npred,umin_mpc, umax_mpc, xlift_min, xlift_max);

% Pre-allocate for storing the data
Xsim = zeros(2*N, Nsim+1);          % Simulated state of controlled FK model
Xsim_unctrl = zeros(2*N, Nsim+1);   % Simulated state of uncontrolled FK model 
Usim = zeros(1, Nsim);              % Controlled Inputs

% Closed-loop simulation
for i = 1:Nsim
    if(mod(i,10) == 0)
        fprintf('Closed-loop simulation, %f %% completed \n', 100 * i / Nsim);
    end

    yref = Xref_idf(:,1);    % reference output without preview (constant reference)
    z = f_BuildKoopmanState(Xsim(:,i));   % create the state of the Koopman linear system via lifting
    
    [u, ~] = kmpc(z,yref); % get the control input
    
    u_disturbance = zeros(N, 1);
    flag_dist = disturbance_idx == i;   % Occurance of the disturbance (and pointer to the array)
    if any(flag_dist)
        disp("Disturbance applied");
        idx_curr_disturb = disturbance_pend_idx(flag_dist);
        u_disturbance(idx_curr_disturb) = disturbance_amp(flag_dist);
    end

    [~,Y_nonlin] = build_and_sim_nonlin_FK_model(fk_params, Xsim(:,i), u(1), u_disturbance, [0 Ts]); % advance in time for one step 
    [~,Y_nonlin_unctrl] = build_and_sim_nonlin_FK_model(fk_params, Xsim_unctrl(:,i), 0, u_disturbance, [0 Ts]); % Advance the uncontrolled system

    % store data
    Xsim(:,i+1) = Y_nonlin(end,:)';
    Usim(i) = u(1);
    Xsim_unctrl(:,i+1) = Y_nonlin_unctrl(end,:)'; 

end

%% Display results
figure;
title('Synchronization to stable equilibrium: Vibration control');
subplot(2,1,1);
hold on;
idx = 1;
for ii = 1:2:2*N
    plot(t, Xsim(ii,:)', 'Linewidth', 1.5, 'Color', color_p(idx));
    idx = idx + 1;
    plot(t, Xsim_unctrl(ii,:)', 'Linewidth', 1, 'Color', color_p(idx), 'LineStyle','--');
end
plot(t, zeros(size(t)),'Linewidth', 3, 'Linestyle', ':');
legend('#1', '#2', '#3','#4','#5', 'Reference');
box on;
grid on;
ylabel('Angle [rad]');
xlabel('Time [s]');
% xlim([0 2]);
ylim([-pi pi]);
subplot(2,1,2);
plot(t(1:end-1), Usim, 'Linewidth', 1.5);
legend('Input');
box on;
grid on;
ylabel('Torque [N.m]');
xlabel('Time [s]');
ylim([-0.15 0.15]);



