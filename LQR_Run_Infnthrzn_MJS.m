% This function implement LQR for MJS with mode switching sequence X, initial state x0, noise w, and
% policy matrices F and W computed with LQR_PolicyCompute_Fnthrzn_MJS.m
% the input u(t)= - W{X(t),t+1} * B{X(t)}' * kron(P(X(t),:), eye(n))*cell2mat(F(:,t+1)) * A{X(t)} * x(:,t);

% A --- Cell type, contains A matrix for each mode
% B --- Cell type, contains B matrix for each mode
% T --- Scalor, Time horizon
% H --- Cost function for states
% L --- Cost function for inputs
% X --- Mode switching sequence, X_{0:T-1}
% x0 --- initial state for mode dynamics
% w --- noise sequence in mode dynamics, w_(0:T-1)
% F, W --- Policy matrices for LQR
% P --- Markov transition matrix

% x --- State sequence x_{0:T}
% u --- Input sequence u_(0:T-1)
% cost --- Cost_{0:T}

function [x, u, cost] = LQR_Run_Infnthrzn_MJS(A, B, T, H, L, X, x0, w, K)
n   = size(A(:,:,1), 1);
m   = size(B(:,:,1), 2);
X0  = X(1);         % Initial mode
X   = X(2:end);     % Mode at time 1 to T-1
w0  = w(:, 1);      % Initial noise
w   = w(:, 2:end);  % noise at time 1 to T-1

u   = zeros(m, T-1);    % u_{1:T-1}
x   = zeros(n, T);      % x_{1:T}, x0 is defined above
cost= zeros(T, 1);      % cost_{1:T}
u0  =  K(:,:,X0)* x0; % Compute u0
cost0 = x0'*H*x0 + u0'*L*u0;
x(:,1)  = A(:,:,X0)*x0 + B(:,:,X0)*u0 + w0; % Compute x1
for t = 1:T-1
    u(:,t)  =  K(:,:,X(t)) * x(:,t);
    cost(t) = x(:,t)'*H*x(:,t) + u(:,t)'*L*u(:,t);
    x(:,t+1)= A(:,:,X(t)) * x(:,t) + B(:,:,X(t)) * u(:,t) + w(:,t);
end
cost(T) = x(:,T)'*H*x(:,T);
x       = [x0, x];    % x_(0:T)
u       = [u0, u];    % u_(0:T-1)
cost    = [cost0; cost]; % cost_{0:T}