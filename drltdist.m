function pMtx = drltdist(alpha)
%   pMtx -- Probability distribution matrix and each row sum to 1
%   alpha -- alpha Parameter matrix, each row serves as the parameter for
%   Dirichelet distribution of corresponding row in pMtx
    pMtx = zeros(size(alpha));
    for i = 1:size(alpha,1)
        for j = 1:size(alpha,2)
            tmp = gamrnd(alpha(i,:),ones(1,size(alpha,2)))';
            pMtx(i,:) = tmp/sum(tmp);
        end
    end
end
