function [x0,z0,x,z] = MJLS_LTIsim(A,B,X,K,T,sigz,sigw)
% Simulate the dynamical system
% x_{k+1} = Ax_k + Bu_k + w_k


% Generate initial dynamics state x0, MC trajectory X_{0:T}, disturbance W_{0:T}
dimX  = size(A, 1);
numModes=size(A, 3);

dimU   = size(B, 2);

X0  = X(1);         % Initial mode
X   = X(2:end);     % Mode at time 1 to T-1
x   = zeros(dimX, T);      % x_{1:T}, x0 is defined above
z  = zeros(dimU, T);      % x_{1:T}, x0 is defined above
x0  =zeros(dimX,1);       % Initial dynamics state distribution
z0 = normrnd(0,sigz,[dimU,1]);
w0 = normrnd(0,sigw,[dimX,1]);

x(:,1)  = (A(:,:,X0)+B(:,:,X0)*K(:,:,X0))*x0...
    +  B(:,:,X0)*z0+w0;
   %

for t = 1:T-1
    z(:,t) = normrnd(0,sigz,[dimU,1]);
    w(:,t)  = normrnd(0,sigw,[dimX,1]);
    x(:,t+1)  = (A(:,:,X(t))+B(:,:,X(t))*K(:,:,X(t)))*x(:,t)...
        + B(:,:,X(t))*z(:,t)+w(:,t);
    %
end

end