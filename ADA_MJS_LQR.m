%% Identification and Adaptive Control of Markov Jump Systems: Sample Complexity and Regret Bounds
%%==========================================================================
%
% Input parameters:
% - A: system matrix; (matrix of dimension dimX x dimX x numModes)
% - B: input matrix; (matrix of dimension dimX x dimU x numModes)
% - Q: state cost matrix; (matrix of dimension dimX x dimX x numModes)
% - R: input cost matrix; (matrix of dimension dimU x dimU x numModes)
% - T: transition matrix; (matrix of dimension numModes x numModes)
%
% Output parameters:
% - Regret 
%
% Last edited: 04 June 2021
%==========================================================================

clear; close all
tic 
dimX_list =  [2];
idxmode_list= [2 ];
simT0=125;
gamma= 2;
numepoch=5;
sig_w_list = [0.01, 0.05, 0.1, 0.12];
for inxdimX=1:numel(dimX_list)
    dimX = dimX_list(inxdimX);  %# of states
    dimU =1;  %# of inputs
    for idxmode=1:numel(idxmode_list)
    numModes = idxmode_list(idxmode);
    numExp = 2;
    eigA=0.3;
    
    %%
    A=zeros(dimX,dimX,numModes);
    B=zeros(dimX,dimU,numModes);
    hA=zeros(dimX,dimX,numModes);
    hB=zeros(dimX,dimU,numModes);
    undeA=zeros(dimX,dimX,numModes);
    undeB=zeros(dimX,dimU,numModes);
    Q = zeros(dimX,dimX,numModes);
    R = zeros(dimU,dimU,numModes);
    hK=zeros(dimU,dimX,numModes); 
    
    for i=1:numModes
        AA= randn(dimX,dimX);                         %%Generate dynamics for each mode
        A(:,:,i) = AA/(abs(eigs(AA,1)))*abs(eigA);    % Scale A so that the leading eigenvalue is eigA
        %A(:,:,i)=randi(2,dimX,dimX);
        %A(:,:,i)= A(:,:,i)./(1.2*max(abs(eig(A(:,:,1)))));
        B(:,:,i) = randn(dimX,dimU);                         % Generate nominal A matrix for each cluster
        tmp = randn(dimX,dimX);
        %Q(:,:,i)  =tmp*tmp';
        Q(:,:,i)  = eye(dimX);
        %tmp = randn(dimU,dimU);
        R(:,:,i)  = eye(dimU);
        %R(:,:,i) = tmp*tmp';
        hK(:,:,i)=randn(dimU,dimX);
    end
    T_e =  drltdist((numModes-1)*eye(numModes, numModes)+1);
    %T_e= drltdist(ones(numModes, numModes));
    
    %% J(hK)
    for idxw = 1:numel(sig_w_list)
        disp(['Start parameter set ', num2str(idxw), '! ', num2str(numel(sig_w_list)), ' in total.'])
        sig_w  = sig_w_list(idxw);  % Used to perturb nominal "A" matrix to generate "A" matrix
        %costt(idxParam,idxParam)=mean((abs(costtt-costtt(1)))/(costtt(1)));
        for idxExp = 1:numExp
            for idxepoch = 1:numepoch
                simT=simT0*gamma^idxepoch;
                sig_z =sig_w/sqrt(simT);
                %%
                t0  = drltdist(ones(1,numModes));    % Initial distribution of MC
                X0  = randsample(1:numModes, 1, true, t0);   % Initial Mode
                tmp = zeros(1, numModes);
                tmp(X0) = 1;
                mc  = dtmc(T_e);
                X  = simulate(mc, simT-1, 'X0', tmp); % X_{0:T-1} % Generate mode switching sequence
                [hA,hB,~]=Inf_MJS_SYID(A,B,hK,simT,X,sig_z,sig_w);
                %%
                hK = cDARE_Inf_MJLS_LQR(hA,hB,Q,R,T_e);
                %%
                w   =  sig_w *randn(dimX,simT);    % w_(0:T-1) Noise in the dynamics
                x0  =1e-4*ones(dimX,1);       % Initial dynamics state distribution
                [~, ~, cost_vec_epoch] = Inf_MJLS_LQR_Cost(hA, hB, simT, Q, R, X, x0, w,hK);
                cost_epoch(idxepoch)=sum(cost_vec_epoch);
                est_exp_cost(inxdimX,idxmode,idxw,idxExp,idxepoch)=sum(cost_epoch);
            end
            %%
            %%
            K = cDARE_Inf_MJLS_LQR(A,B,Q,R,T_e);
            t0  = drltdist(ones(1,numModes));    % Initial distribution of MC
            X0  = randsample(1:numModes, 1, true, t0);   % Initial Mode
            tmp = zeros(1, numModes);
            tmp(X0) = 1;
            mc  = dtmc(T_e);
            X  = simulate(mc, simT-1, 'X0', tmp); % X_{0:T-1} % Generate mode switching sequence
            x0  = zeros(dimX,1);       % Initial dynamics state distribution
            w   = sqrt(sig_w)*randn(dimX,simT);    % w_(0:T-1) Noise in the dynamics
            [~, ~, exa_exp_mode_cost] = Inf_MJLS_LQR_Cost(A, B, simT, Q, R, X, x0, w,K);
            exa_exp_cost(inxdimX,idxmode,idxw,idxExp)=sum(exa_exp_mode_cost);
        end
        regret= est_exp_cost-exa_exp_cost;        
    end
    end
end
toc

%filename = 'regret_';
%save(['regret/',filename, datestr(now, "yyyymmddHHMMSS"), '.mat']);



Marksize  = 12;
linewidth = 6;
Fontsize  = 50;
Err1=squeeze(regret(1,1,1,:,:));
Err2=squeeze(regret(1,1,2,:,:));
Err3=squeeze(regret(1,1,3,:,:));
Err4=squeeze(regret(1,1,4,:,:));
simT_list  = [250 500 1000 2000 4000];

% line
figure;
     shadedErrorBar(simT_list, Err1, {@mean,@std}, 'lineprops', '-m'); hold on
    % shadedErrorBar(simT_list, Err2, {@mean,@std}, 'lineprops', '-b');
     shadedErrorBar(simT_list, Err2, {@mean,@std}, 'lineprops', '--g');
     shadedErrorBar(simT_list, Err3, {@mean,@std}, 'lineprops', ':c');
     shadedErrorBar(simT_list, Err4, {@mean,@std}, 'lineprops', '-.k');
h1 = semilogy(simT_list,mean(Err1), 'marker','o','markersize',Marksize,'markerfacecolor','m');
     semilogy(simT_list,mean(Err1),'-m','linewidth',linewidth,'marker','o','markersize',Marksize,'markerfacecolor','m');
%h2 = semilogy(simT_list,mean(Err2),'marker','d','markersize',Marksize,'markerfacecolor','b');
%     semilogy(simT_list,mean(Err2),'b','linewidth',linewidth,'marker','d','markersize',Marksize,'markerfacecolor','b');
h3 = semilogy(simT_list,mean(Err2),'marker','s','markersize',Marksize,'markerfacecolor','g');
     semilogy(simT_list,mean(Err2),'--g','linewidth',linewidth,'marker','s','markersize',Marksize,'markerfacecolor','g');
h4 = semilogy(simT_list,mean(Err3),'marker','p','markersize',Marksize,'markerfacecolor','c');
     semilogy(simT_list,mean(Err3),':c','linewidth',linewidth,'marker','p','markersize',Marksize,'markerfacecolor','c');
h5 = semilogy(simT_list,mean(Err4),'marker','>','markersize',Marksize,'markerfacecolor','k');
     semilogy(simT_list,mean(Err4),'-.k','linewidth',linewidth,'marker','>','markersize',Marksize,'markerfacecolor','k');


% limits  
set(gca,'FontSize',30)
grid on
grid minor
%ylim([10,1e+3])
xlim([250, 4000])
xticks([  500 1000 2000 4000])
xticklabels({ '500', '1000', '2000', '4000'})
ytickformat('%1.0e');
%set(gca,'TickLabelInterpreter','latex','FontSize',Fsize);
     
% label
xlabel('$T$','Interpreter','latex','FontSize',Fontsize,'FontWeight','bold','Color','k');
ylabel('Regret','Interpreter',...
      'latex','FontSize',Fontsize,'FontWeight','bold','Color','k');

% figure size
set(gcf,'Position',[250 150 300 320]);
%legend('$n = 5$','$n = 10$', '$n =15$', '$n =20$',...
   %    'Interpreter','latex', 'Location','Northeast','FontSize',Fontsize);

legend('$\sigma_w = 0.01$','$\sigma_w = 0.05$', '$\sigma_w = 0.1$', '$\sigma_w = 0.2$',...
       'Interpreter','latex', 'Location','Northeast','FontSize',Fontsize);