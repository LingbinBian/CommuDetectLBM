function [connectpara_infer]=element_modelling(adj_cell,N,S)
% This function infers the block mean and variance of communities for a single subject.
%
% Input:
%    adj_cell: cell of connectivity
%    esti_label: estimation of labels
%    K: number of communities
%    N: number of nodes
%    S: replication number
% Output:
%    localpara_infer: a struct containing infered parameters of the network
%
% Version 1.0
% 11-April-2021
% Copyright (c) 2021, Lingbin Bian
% -------------------------------------------------------------------------

% hyper-parameters of prior
nu=3;
rho=0.02;
xi=0;
kappa_sq=1;
S_subj=100;     % number of subjects
 
Preci=cell(S,1); % connectivity precision matrix
Vari=cell(S,1); % connectivity variance
Mean=cell(S,1); % connectivity mean

for s=1:S
   % sampling model parameters from the posterior conditional on z 
   [Preci{s},Vari{s},Mean{s}]=sampling_para(adj_cell,S_subj); 
end

connectpara_infer=struct('Preci', Preci,...
                 'Vari',Vari,...
                 'Mean',Mean,...
                 'S',S);
                          
% Nested functions---------------------------------------------------------
%--------------------------------------------------------------------------

% Sampling mean and variance from the posterior p(u,sigma^2|x,z)
    function [precision,vari,mean] = sampling_para(x,S_subj)  %O(K*N,K^2,N^2)
         
        w_sum=zeros(N,N);    % sum of weights
        w_sumsq=zeros(N,N);  % sum of square

        for i=1:N
            for j=1:N
                    w_sum(i,j)=sum(x{i,j});
                    w_sumsq(i,j)=sumsqr(x{i,j});           
            end
        end
        vari=zeros(N,N);
        mean=zeros(N,N);
        precision=zeros(N,N);
        for i=1:N
            for j=1:N
                nu_s=nu+S_subj;
                rho_s=(xi^2)/(kappa_sq)+w_sumsq(i,j)+rho-((xi+w_sum(i,j)*kappa_sq)^2)/(1/kappa_sq+S_subj);
                precision(i,j)=gamrnd(0.5*nu_s,1/(0.5*(rho_s)));
                vari(i,j)=1/precision(i,j);    % inverse gamma
            end
        end
        vari(isnan(vari))=0;
        for i=1:N
            for j=1:N              
                mean(i,j)=normrnd((xi+w_sum(i,j)*kappa_sq)/(1+S_subj*kappa_sq),sqrt((kappa_sq/(1+S_subj*kappa_sq))*vari(i,j)));
            end
        end  
    end
%--------------------------------------------------------------------------
end

