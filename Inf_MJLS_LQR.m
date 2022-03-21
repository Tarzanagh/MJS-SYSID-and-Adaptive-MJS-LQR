function [L,Pl,Exx] = Inf_MJLS_LQR(A,B,H,Q,R,T)
%==========================================================================
% Inf_MJLS_LQR: computes the regulator gain for static feedback
%               control of Markov jump linear systems.
%
% Input parameters:
% - A: system matrix; (matrix of dimension dimX x dimX x numModes)
% - B: input matrix; (matrix of dimension dimX x dimU x numModes)
% - H: noise matrix; (matrix of dimension dimX x dimW x numModes)
% - Q: state cost matrix; (matrix of dimension dimX x dimX x numModes)
% - R: input cost matrix; (matrix of dimension dimU x dimU x numModes)
% - T: transition matrix; (matrix of dimension numModes x numModes)
%
% Output parameters:
% - L: regulator gain; (matrix of dimension dimU x dimX)
%
%
% last edited: 30 Jan 2021
%==========================================================================
itermax=100;
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
stepsize=0.25;
Exx = zeros(dimX,dimX,numModes);
La = zeros(dimX,dimX,numModes);
L = zeros(dimU,dimX,numModes);
for i = 1:numModes
    Exx(:,:,i) =1e-4*eye(dimX,dimX);
    La(:,:,i) = 1e-4*eye(dimX,dimX);
end
oldL = L;
% END: initialize

% iteration until convergence
for k = 1:itermax
    P = OpEpsilon(La,T);
    
        % compute L
    for i = 1:numModes
        XRBPB = kron(Exx(:,:,i),R(:,:,i)+B(:,:,i)'*P(:,:,i)*B(:,:,i));
        BPAX =  B(:,:,i)'*P(:,:,i)*A(:,:,i)*Exx(:,:,i);
        L(:,:,i) =reshape(-XRBPB\BPAX(:),[dimU,dimX]);
    end
    % END: compute L
    
%      compute L via gradient descent
%       for i = 1:numModes
%          XRBPB= kron(Exx(:,:,i),R(:,:,i)+B(:,:,i)'*P(:,:,i)*B(:,:,i));
%          GRAD_L = 2*((R(:,:,i)+B(:,:,i)'*P(:,:,i)*B(:,:,i))*oldL(:,:,i)...
%                 + B(:,:,i)'*P(:,:,i)*A(:,:,i))*Exx(:,:,i);
%          L(:,:,i) =oldL(:,:,i)-stepsize*reshape(XRBPB\GRAD_L(:),[dimU,dimX]);
%       end
%       END: compute L
    
    Exx_ = RecursionExx(Exx,L,A,B,H,T,mode);
    La_ = RecursionLambda(La,L,A,B,Q,R,T,mode);
    
    % check convergence
    if(k>1)
        if(sum(sum(sum(abs(oldL- L)))) < tol)
            break
        end
    end
    % END: check convergence
    
    Exx = Exx_;
    La = La_;
    oldL = L;
end

Pl = zeros(dimX,dimX,numModes);
for i = 1:numModes
    Pl(:,:,i) =A(:,:,i)'*P(:,:,i)*A(:,:,i)+Q(:,:,i)...
           - (A(:,:,i)'*P(:,:,i)*B(:,:,i))*L(:,:,i);
end
% END: % iteration until convergence
end % function MJLS_StaticOutputComputeRegulatorGain

%% RecursionExx
function outExx = RecursionExx(Exx,L,A,B,H,T,mode)
outExx = zeros(size(Exx));
numModes = numel(mode);

for j = 1:numModes
    for i = 1:numModes
        outExx(:,:,j) = outExx(:,:,j) +...
            T(i,j)*((A(:,:,i)+ B(:,:,i)* L(:,:,i))*Exx(:,:,i)*(A(:,:,i)...
            + B(:,:,i)* L(:,:,i))' + mode(i)*H(:,:,i)*H(:,:,i)');
    end
    outExx(:,:,j) = (outExx(:,:,j) + outExx(:,:,j)')/2;
end
end % function RecursionExx

%% RecursionLambda
function outP = RecursionLambda(P,L,A,B,Q,R,T,mode)
outP = zeros(size(P));
numModes = numel(mode);

eP = OpEpsilon(P,T);

for j = 1:numModes
    outP(:,:,j) = (A(:,:,j)...
        + B(:,:,j)*L(:,:,j))'*eP(:,:,j)*(A(:,:,j)...
        + B(:,:,j)*L(:,:,j)) + Q(:,:,j) + L(:,:,j)'*R(:,:,j)*L(:,:,j);
    outP(:,:,j) = (outP(:,:,j) + outP(:,:,j)')/2;
end
end % function RecursionLambda

%% OpEpsilon
function out = OpEpsilon(M,T)
[numR,numC,nModes] = size(M);
out = zeros(numR,numC,size(T,1));
for i = 1:nModes
    for j = 1:nModes
        out(:,:,i) = out(:,:,i) + T(i,j)*M(:,:,j);
    end
end
end % function OpEpsilon