function [U, X, MPC_ctr] = osqp_MPC_controller(A,B,C,Q,R,QN,N,u_min, u_max, z_min, z_max)
% MPC in dense formulation
% Authors: Milan Korda, Hassan Arbabi, Loi Do
% doloi@fel.cvut.cz
% @2022

% OUTPUT:
% Tracking MPC controller MPC_ctr(x0,yr)

% Dynamics:
% x^+ = A*x + B*u
% y   = C*x

% Cost:
% J =  1/2*(y_N - yr_N)'*Q_N*(y_N - yr_N) + 1/2*sum_{i=0:N-1} [(y_i - yr_i)'*Q*(y_i - yr_i) + u_i'*R*u_i]

% INPUTS:
% Xlb, Xub:  n x 1 or n x N matrix
% Ulb, Uub: m x 1 vector or m x N matrix

% Example: [~, ~, mpcCont ]  = osqp_MPC_controller(A,B,C,d,Q,R,QN,N,u_min, u_max, z_min, z_max)
% To get just the controller
% MPC_ctr(x0,yr) generates the control input

n = size(A,1); % Number of states
m = size(B,2); % Number of control inputs
p = size(C,1); % Number of outputs

% ------------------------------
[Ab_bar, Ab, Bb_bar, Bb] = stackedMPCmatrices(A,B,N);
Qb = sparse(blkdiag(kron(eye(N-1), Q), QN));    % Size: Np*output_size x Np*output_size
Rb = sparse(kron(eye(N), R));                   % Size: Np*input_size x Np*input_size
Cb = sparse(kron(eye(N), C));                   % Size: Np*output_size x Np*output_size
% ------------------------------

% % Bounds on the states and inputs
% Stack for all prediction steps
z_min = repmat(z_min, N+1, 1); % +1 for z_0
z_max = repmat(z_max, N+1, 1); % +1 for z_0
u_max = repmat(u_max, N, 1);
u_min = repmat(u_min, N, 1);

b_min = [z_min; u_min];
b_max = [z_max; u_max];

% % ------------ Gather all matrices and vectors for warm start of the optimization program ------------
%      min      1/2*x'Hx + f'x
%      s.t.     l <= A_ineq*x <= u

H = Bb'*Cb'*Qb*Cb*Bb + Rb;

y_ref_curr = zeros(N*p,1);
x0 = zeros(n,1);
z_curr = x0;

H = sparse((H+H')/2); % symetrize (in case of numerical errors)
A_ineq =  sparse(Bb_bar);

l = b_min - Ab_bar*z_curr;
u = b_max - Ab_bar*z_curr;

% Precomputed matrices to compute linear part of the cost 
M_z = Ab'*(Cb'*Qb*Cb)*Bb; % Multiplying lifted state
M_ref = -Qb*Cb*Bb; % Multiplying output reference
f = z_curr'*M_z + y_ref_curr'*M_ref;

% Build the controller
disp('Building controller with osqp');
QP_solver = osqp;
QP_solver.setup(H, f, A_ineq, l, u, 'warm_start', true,  'verbose', false, 'eps_abs', 1e-04, 'eps_rel', 1e-04);
% QP_solver.setup(H, f, A_ineq, l, u, 'warm_start', true,  'verbose', false, 'polish', 1);
solution = QP_solver.solve();

MPC_ctr = @(z_curr,y_ref_curr)(mpcController_osqp(z_curr,y_ref_curr,QP_solver,N,Ab_bar,A_ineq,b_min,b_max,M_z,M_ref, C, Q));   

U = solution.x;
X = [x0;Ab*x0+Bb*U];

end




function [U, optval] = mpcController_osqp(z_curr, y_ref_curr,QP_solver,N,Ab_bar,A_ineq,b_min,b_max,M_z,M_ref, C, Q)
% INPUTS: 
% - z_curr: lifted state at current time
% - y_ref_curr: reference signal (on output of the lifted system) with preview 

% PARAMETERS
% - N: prediction horizon
% - A_ineq: defines constraints:    l <= A_ineq*x <= b
% - b_min: [z_min_stacked, u_min_stacked];      max/min bounds on lifted states and input stacked for all prediction steps
% - b_max: [z_max_stacked, u_max_stacked];       
% - M_z, M_ref: precomputed matrices to compute linear part of the cost

% M_z = (2*Ab'*Cb'*Qb*Cb*Bb)';
% M_ref = (-2*Qb*Cb*Bb)';

if (size(y_ref_curr,2) == 1)
    % Without preview, stack one-step reference:
    y_ref_curr_stacked = repmat(y_ref_curr, N, 1);
else
    % With preview, reshape into vector
    y_ref_curr_stacked = y_ref_curr(:);
end

%      min      1/2*x'Hx + f'x
%      s.t.     l <= A_ineq*x <= b

% Linear part of the objective function
f = z_curr'*M_z + y_ref_curr_stacked'*M_ref;

% ---- Solve QP ----
% Update initial state
l = b_min - Ab_bar*z_curr;
u = b_max - Ab_bar*z_curr;
QP_solver.update('q', f, 'l', l, 'u', u);

% Check solver status
solution = QP_solver.solve();
if ~strcmp(solution.info.status, 'solved')
    warning('OSQP did not solve the problem!');
end

% Give solution
U = solution.x;

e = (y_ref_curr(:,1) - C*z_curr) ;
optval = solution.info.obj_val + e'*Q*e;


end






