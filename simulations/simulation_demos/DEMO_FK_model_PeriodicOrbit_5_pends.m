%% DEMO: synchronization to the Periodic orbit
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

Qy = kron(eye(N), diag([10,0.5])); % Penalization of system's outputs
Qy_spring = kron(L, [fk_params.k, 0; 0, 0]);
Qy_dissip = kron(L, [0, 0; 0, fk_params.b]);
Qy = Qy + 3000*(Qy_dissip + Qy_spring);

R = 1;

Npred = 50; % Prediction horizon for KMPC

nu = size(R,2); % input dimension
umin_mpc = umin.*ones(nu,1);
umax_mpc = umax.*ones(nu,1);

%%% Identification
Npred_idf = 1;
Ntraj_idf = 100;    % Number of trajectories for identification
Tsim_idf = 5;    % Simulation time of a single trajectory for identification  
X0_idf =  repmat([0;0],N,Ntraj_idf); % Initial states for all trajectories starting at origin
t_idf = (0:Ts:(Tsim_idf + (Npred_idf+10)*Ts))'; % Time for reference trajectory
v_ref = 15;     % Angular speed of the reference trajectory
Xref_idf = repmat([t_idf*v_ref, ones(size(t_idf))*v_ref]', N, 1); % Reference for all pendulums and for whole time interval of identification

Qy_idf = kron(eye(N), diag([100,0.1]));
R_idf = 1;

u_MPC_random_sigma = 0.5;  % Perturbation to the input used for identification

x_min = -10000*ones(2*N,1); % Lower bounds on states
x_max = 10000*ones(2*N,1); % Upper bounds on states

%%% Final simulation
n = 6*N;  % lifted state dimension
xlift_min = -10000*ones(n,1); % Lower bounds on lifted states
xlift_max = 10000*ones(n,1); % Upper bounds on lifted states

Tsim = 10; % Simulation length of the final simulation
Nsim = Tsim/Ts;
t = (0:Ts:(Tsim+(Npred+10)*Ts))'; % prolong the reference signal for KMPC preview
Xref = repmat([t*v_ref, ones(size(t))*v_ref]', N, 1);

%% ************************PART I: Identification ********************** %%
% Create controller based on linearized FK model.
[sys_d, ~] = create_linearized_FK_model(fk_params, Ts); % Linearized system around system's origin 

% Build MPC controller with linearized system
[~,~,kmpc]  = osqp_MPC_controller(sys_d.A,sys_d.B,sys_d.C,Qy_idf,R_idf,Qy_idf,Npred_idf,umin_mpc, umax_mpc, x_min, x_max);
f_kmpc = @(x, r) kmpc(x,r) + normrnd(0, u_MPC_random_sigma);    % Add random perturbation

% Collecting closed-loop data for EDMD
f_collect_data = 1; % Flag: 1: collect data, 0: load saved data
if f_collect_data == 1
    [X,Y,U] = CollectData_FK_model_closed_loop(Ntraj_idf, Tsim_idf, Ts, f_kmpc, Xref_idf, X0_idf, fk_params);
    f_name = ['./simulations/data/DEMO_data_PeriodicTraj_5_pends_numTraj-' , num2str(Ntraj_idf), '.mat'];
    save(f_name, 'X', 'Y', 'U');
else
    load('./simulations/data/DEMO_data_PeriodicTraj_5_pends_numTraj-100.mat');  
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
Usim = zeros(1, Nsim);              % Controlled Inputs

% Closed-loop simulation
for i = 1:Nsim
    if(mod(i,10) == 0)
        fprintf('Closed-loop simulation, %f %% completed \n', 100 * i / Nsim);
    end

    yref = Xref(:,i:i+Npred-1);    % Current reference with preview
    z = f_BuildKoopmanState(Xsim(:,i));   % create the state of the Koopman linear system via lifting
    
    [u, ~] = kmpc(z,yref); % get the control input
    u_disturbance = zeros(N, 1); % No disturbance

    [~,Y_nonlin] = build_and_sim_nonlin_FK_model(fk_params, Xsim(:,i), u(1), u_disturbance, [0 Ts]); % advance in time for one step 

    % store data
    Xsim(:,i+1) = Y_nonlin(end,:)';
    Usim(i) = u(1);

end

%% Display results
figure;
title('Synchronization to Periodic orbit');


subplot(2,2,1);
hold on;
idx = 1;
for ii = 1:2:2*N
    plot(t(1:size(Xsim, 2)), Xsim(ii,:)', 'Linewidth', 1.5, 'Color', color_p(idx));
    idx = idx + 1;
end
plot(t, Xref(1,1:numel(t)),'Linewidth', 3, 'Linestyle', ':');
box on;
grid on;
ylabel('Angle [rad]');
xlabel('Time [s]');


subplot(2,2,2);
hold on;
idx = 1;
for ii = 2:2:2*N
    plot(t(1:size(Xsim, 2)), Xsim(ii,:)', 'Linewidth', 1.5, 'Color', color_p(idx));
    idx = idx + 1;
end
plot(t, Xref(2,1:numel(t)),'Linewidth', 3, 'Linestyle', ':');
box on;
grid on;
ylabel('Angle [rad]');
xlabel('Time [s]');



% Input
subplot(2,2,3:4);
plot(t(1:numel(Usim)), Usim, 'Linewidth', 1.5);
legend('Input');
box on;
grid on;
ylabel('Torque [N.m]');
xlabel('Time [s]');


