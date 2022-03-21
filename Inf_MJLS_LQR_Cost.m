function [x, u, cost] = Inf_MJLS_LQR_Cost(A, B, T, Q, R, X, x0, w,K)
       
 % Generate initial dynamics state x0, MC trajectory X_{0:T}, disturbance W_{0:T}
n   = size(A, 1);
m   = size(B, 2);
X0  = X(1);         % Initial mode
X   = X(2:end);     % Mode at time 1 to T-1
w0  = w(:, 1);      % Initial noise
w   = w(:, 2:end);  % noise at time 1 to T-1

u   = zeros(m, T-1);    % u_{1:T-1}
x   = zeros(n, T);      % x_{1:T}, x0 is defined above
cost= zeros(T, 1);      % cost_{1:T}
u0  =  K(:,:,X0) * x0; % Compute u0
cost0 = x0'*Q(:,:,X0)*x0 + u0'*R(:,:,X0)*u0;
x(:,1)  = A(:,:,X0)*x0 + B(:,:,X0)*u0 + w0; % Compute x1
for t = 1:T-1
    u(:,t)  =  K(:,:,X0) * x(:,t);
    cost(t) = x(:,t)'*Q(:,:,X(t))*x(:,t) + u(:,t)'*R(:,:,X(t))*u(:,t);
    x(:,t+1)= A(:,:,X(t)) * x(:,t) + B(:,:,X(t)) * u(:,t) + w(:,t);
end
cost(T) = x(:,T)'*Q(:,:,X0)*x(:,T);
x       = [x0, x];    % x_(0:T)
u       = [u0, u];    % u_(0:T-1)
cost    = [cost0; cost]; % cost_{0:T}