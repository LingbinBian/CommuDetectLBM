% This script tests the MCMC sampler (Gibbs & M3)
%
% Version 1.0
% 15-March-2020
% Copyright (c) 2020, Lingbin Bian
% -------------------------------------------------------------------------
clear
clc
% -------------------------------------------------------------------------
% Load data 
% Data type
datatype=1;   % 1: real data, 0: synthetic data

if datatype==1
    subjects=load('subject.txt');
    localmin_t=[41,76,107,140,175,207,238,278,306,333,375];
    K_min=[6,6,6,6,6,6,6,6,6,6,6];  % make an assumption of K
    W=10;
elseif datatype==0
    subjects=load('synthetic_id.txt');
    localmin_t=[36,66,92,116,146];
    K_min=[4,5,3,5,4];
    W=10;
end

N_subj=100;
N=35;
S=200;
L_localmin=length(localmin_t);

group_adj=cell(N_subj,1); % group of adjacency matrix

esti_grouplabel=zeros(N,L_localmin);
ave_adj=cell(1,L_localmin);  % averaged adjacency matrix

% Local averaged adjacency matrix
for s=1:N_subj
    fprintf('Adjacency of subject: %d\n',s)
    subid=num2str(subjects(s));
    [group_adj{s,1},true_latent]=local_adj(datatype,subid,1,0.5623,localmin_t,K_min,W);
end

adj_mean=zeros(N_subj,1);
for t=1:L_localmin
    for i=1:N
        for j=1:N
            for s=1:N_subj
                adj_mean(s,1)=group_adj{s,1}{1,t}(i,j);
            end
            ave_adj{1,t}(i,j)=mean(adj_mean);
        end
    end
end
% -------------------------------------------------------------------------
% Adjacency matrix
K=5;
W=10; % half of window size
M=10; % margin size
N=35; 

x=ave_adj{1,4};   % single adjacency matrix


%--------------------------------------------------------

% -------------------------------------------------------------------------
% Prior hyper-parameters 
Itera=1000;
nu=3;
rho=0.02;
xi=0;
kappa_sq=1;

priorpra=struct('nu',nu,...
                'rho',rho,...
                'xi',xi,...
                'kappa_sq',kappa_sq,...
                'Itera',Itera); 
            
% -------------------------------------------------------------------------
% MCMC Allocation sampler

% Multiple simulations
N_sim=2;   % number of simulations
chain_simul=zeros(Itera,1,N_sim);
tic
for i=1:N_sim
 [Z_chain,prob_chain,accep_r_array,K_chain]=MCMC_allocation(x,K,priorpra);
 chain_simul(:,:,i)=prob_chain';
end
% toc
% figure 
% scatter(1:Itera,accep_r_array(1:Itera),10,'filled')
% title('Acceptance ratio','fontsize',16) 
% xlabel('Iteration','fontsize',14)
% ylabel('r','fontsize',14)
% 
% % -----------------------------------------------------------------------
% % Convergence diagnostic
% PSRF=zeros(1,Itera);  % Potential scale reduction factor
% for k=2:Itera   
%     [PSRF(1,k),NEFF(1,k),V(1,k),W(1,k),B(1,k)]=psrf(chain_simul(1:k,:,:)); 
% end
% figure
% plot(20:Itera,PSRF(20:Itera),'Linewidth',2.0)
% title('Potential Scale Reduction Factor of p(z)','fontsize',16) 
% xlabel('Iteration','fontsize',14)
% ylabel('PSRF','fontsize',14)


% Test esti_community
[Z_chain,K_chain,z,K]=esti_community(x)

% -------------------------------------------------------------------------
% plot adjacency matrix
figure
imagesc(x)
title('Adjacency matrix','fontsize',16) 
xlabel('Node number','fontsize',14)
ylabel('Node number','fontsize',14)
colormap(hot);
colorbar;
set(gca,'position',[0.1,0.2,0.7,0.6] );
set(gca,'fontsize',14)
set(gcf,'unit','normalized','position',[0.1,0.2,0.16,0.27]);
set(gca, 'linewidth', 1.2, 'fontsize', 16, 'fontname', 'times')
%--------------------------------------------------------------------------
% plot Markov chain
figure
imagesc(Z_chain)
title('MCMC allocation sampler','fontsize',16) 
xlabel('Iteration','fontsize',14)
ylabel('Latent labels','fontsize',14)
%colormap(parula(max(max(Z_chain))));
%colorbar_community(max(max(Z_chain)));
colormap(parula(max(K_chain)));
colorbar_community(max(K_chain));
set(gca,'position',[0.07,0.2,0.8,0.65] );
set(gca,'fontsize',14)
set(gcf,'unit','normalized','position',[0.3,0.2,0.4,0.27]);
set(gca, 'linewidth', 1.2, 'fontsize', 16, 'fontname', 'times')

figure
scatter(1:Itera,K_chain(1,:),5,'filled','d')
title('Markov chain of K','fontsize',16) 
xlabel('Iteration','fontsize',14)
ylabel('K','fontsize',14)
set(gca,'position',[0.07,0.2,0.8,0.65] );
set(gca,'fontsize',14)
set(gcf,'unit','normalized','position',[0.3,0.2,0.4,0.27]);
set(gca, 'linewidth', 1.2, 'fontsize', 16, 'fontname', 'times')
% cmap = jet(3)
% %cmap = flipud(cmap(1:4,:));
% %cmap(1,:) = [1,1,1];
% colormap(cmap);
% colorbar
%--------------------------------------------------------------------------

