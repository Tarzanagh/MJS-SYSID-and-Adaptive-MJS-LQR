function [K,outP] = cDARE_Inf_MJLS_LQR(A,B,Q,R,T)
%==========================================================================
% Inf_MJLS_LQR: computes the feedback matrices for the MJS-LQR problem .
%
% Input parameters:
% - A: system matrix; (matrix of dimension dimX x dimX x numModes)
% - B: input matrix; (matrix of dimension dimX x dimU x numModes)
% - Q: state cost matrix; (matrix of dimension dimX x dimX x numModes)
% - R: input cost matrix; (matrix of dimension dimU x dimU x numModes)
% - T: transition matrix; (matrix of dimension numModes x numModes)
%
% Output parameters:
% - K: regulator gain; (matrix of dimension dimU x dimX)
%
%
% LastUpdate: 25 Feb 2021
%==========================================================================
itermax=1e4;
tol=1e-6;

[dimX,dimU,numModes] = size(B);
mode = zeros(numModes,1);
mode(1) = 1;

for k = 1:itermax
    mode_ = mode;
    mode = T'*mode_;
    if(sum(abs(mode-mode_)) < tol)
        break
    end
end
if(sum(abs(mode-mode_)) > tol)
    error('Mode did not converge')
end

% initialization
for i = 1:numModes
    outP(:,:,i) = Q(:,:,i);
end
oldP = OpEpsilon(outP,T);
% END: initialize
% iteration until convergence
for k = 1:itermax
    numModes = numel(mode);
    for i = 1:numModes
        RBPB(:,:,i) = R(:,:,i)+B(:,:,i)'*oldP(:,:,i)*B(:,:,i);
        outP(:,:,i) = A(:,:,i)'*oldP(:,:,i)*A(:,:,i)+Q(:,:,i)...
            - (A(:,:,i)'*oldP(:,:,i)*B(:,:,i))...
            * (RBPB(:,:,i)\(B(:,:,i)'*oldP(:,:,i)*A(:,:,i))); % cDARE
    end
    % END: compute P
    % check convergence
    P=OpEpsilon(outP,T);
    if(k>1)
        if sum(sum(sum(abs(oldP- P))))  < tol
            break
        end
    end
    % END: check convergence
    oldP = P;
end

% compute K
K = zeros(dimU,dimX,numModes);
for i = 1:numModes
    RBPB(:,:,i) = R(:,:,i)+B(:,:,i)'*P(:,:,i)*B(:,:,i);
    K(:,:,i) = -RBPB(:,:,i)\(B(:,:,i)'*P(:,:,i)*A(:,:,i));
end
% END: % iteration until convergence
end
%
%% OpEpsilon
function out = OpEpsilon(M,T)
[numR,numC,nModes] = size(M);
out = zeros(numR,numC,size(T,1));
for i = 1:nModes
    for j = 1:nModes
        out(:,:,i) = out(:,:,i) + T(i,j)*M(:,:,j);
    end
end
end