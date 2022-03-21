function X =Gene_Mode_Swit(T_e,simT)     
numModes= size(T_e,2);
for ii = 1:numModes
            t0  = drltdist(ones(1,numModes));    % Initial distribution of MC
            X0  = randsample(1:numModes, 1, true, t0);   % Initial Mode
            tmp = zeros(1, numModes);
            tmp(X0) = 1;
            mc  = dtmc(T_e);
            X  = simulate(mc, simT-1, 'X0', tmp); % X_{0:T-1} % Generate mode switching sequence
end
end