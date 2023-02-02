function [A,B,C,BuildKoopmanState]= SystemID_via_EDMD_FK(X,Y,U)
% Extended library with: sin(x), cos(x), v*sin(x), v*cos(x)


%% Extended Dynamic Mode Decomposition

% Add sine of positions
sine_X = @(y) sin(y(1:2:end, :));

% Add cosine of positions
cosine_X = @(y) cos(y(1:2:end, :));

% Add v*sin(x) and v*cos(x)
v_sine_X = @(y) y(2:2:end, :).*sine_X(y);
v_cosine_X = @(y) y(2:2:end, :).*cosine_X(y);

% Stack observables
Xp = [X; sine_X(X); cosine_X(X); v_sine_X(X); v_cosine_X(X)];
Yp = [Y; sine_X(Y); cosine_X(Y); v_sine_X(Y); v_cosine_X(Y)];
Up = U;

% compute the Koopman linear system 
tic
disp('Running EDMD ...')
W = [Xp;Up]*[Xp;Up]';
V = Yp*[Xp;Up]';
M = V*pinv(W);

Nlift = size(Xp,1);
A = M(:,1:Nlift);
B = M(:,Nlift+1:end);

nm= size(X,1);
C = zeros(nm,size(A,1));
C(:,1:nm)=eye(nm);
toc

%% Build Koopman State from the current state 'x' for Closed-loop simulations
BuildKoopmanState = @(x) [x; 
    sine_X(x); cosine_X(x);
    v_sine_X(x); v_cosine_X(x)];

end
