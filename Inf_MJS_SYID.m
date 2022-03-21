function [hA,hB,hG,X]=Inf_MJS_SYID(A,B,K,T,X,sigz,sigw)
dimX  = size(A, 1);
numModes=size(A, 3);
dimU   = size(B, 2);
hA=zeros(dimX,dimX,numModes);
hB=zeros(dimX,dimU,numModes);

%%MJS
[x0,z0,x,z] = MJLS_LTIsim(A,B,X,K,T,sigz,sigw);

H=[];
Y=[];
for i=1:numModes
    if X(1)==i
        H=[ x0'   z0'] ;
        Y=[ x(:,1)' ] ;
    end
    for t=1:T-1
        if X(t+1)==i
            if  isempty(H)
                H=[ x(:,t)'    z(:,t)' ] ;
                Y=[ x(:,t+1)' ] ;
            else
                H=[H;  x(:,t)'    z(:,t)' ] ;
                Y=[Y; x(:,t+1)' ] ;
            end
        end
    end
    hG(:,:,i) = pinv(H)*Y;
    hGi =hG(:,:,i);
    hB(:,:,i)=hGi(dimX+1:end,:)';
    hA(:,:,i)=hGi(1:dimX,1:dimX)-hB(:,:,i)*K(:,:,i);
    H=[];
    Y=[];

end