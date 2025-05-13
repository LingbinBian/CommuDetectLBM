function [z_esti,K,z_chain,K_chain] = CommuDetectLBM(x)
% Community detection for a single adjacency matrix
% Latent block model (LBM) with Markov chain Monte Carlo (MCMC) allocation sampler
%
% Input:
%         x: NxN adjacency matrix 
%
% Output:
%         z_esti: a vector of community memberships
%         K: number of communities
%         z_chain: a generated Markov chain updating community memberships (N x Itera matrix)
%         K_chain: a generated chain of the number of communities
%
% In this function, the users can set several parameters of the
% model. These parameters are divided into two sets. The first sets of
% parameters are the prior hyper-parameters of the latent block model
% including nu,rho,xi, and kappa_sq, where xi and kappa_sq are the parameters of the distribution of block mean
% and nu and rho are the parameters of the distribution of block variance.
% The second sets of parameters are related to the generated Markov chain,
% K_ini is the initial assumption of the number of communities at the first
% step of the chain, Itera is the total number of steps of the chain,
% bunin_ite is the burn-in steps, S is the thinning value and
% autocorre_time is the autocorrelation.

% Parameter setting:
%--------------------------------------------------------------------------
% example (default)
% Prior hyper-parameters of latent block model
% nu=3;
% rho=0.02;
% xi=0;
% kappa_sq=1;

% parameters of Markov chain
% K_ini=5;
% Itera=2000;
% burnin_ite=500; 
% S=400; % sample size
% autocorre_time=3;
%--------------------------------------------------------------------------
% Version 1.0
% 10-April-2025
% Copyright (c) 2025, Lingbin Bian, 
% Contact via email: lingbin.bian@gmail.com
%--------------------------------------------------------------------------
% Parameter setting
%--------------------------------------------------------------------------
% Prior hyper-parameters of latent block model
nu=3;
rho=0.02;
xi=0;
kappa_sq=1;

% parameters of Markov chain
K_ini=5;
Itera=2000;
burnin_ite=500;
S=400;
autocorre_time=3;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

priorpra=struct('nu',nu,...
                'rho',rho,...
                'xi',xi,...
                'kappa_sq',kappa_sq,...
                'Itera',Itera); 
            
N=length(x);
% MCMC allocation sampler
[z_chain,prob_chain,accep_r_array,K_chain]=MCMC_allocation_sampler(x,K_ini,priorpra);


K_pos=zeros(1,S);
z_pos=zeros(N,S);

for s=1:S
   z_pos(:,s)=z_chain(:,burnin_ite+autocorre_time*s);
   K_pos(1,s)=K_chain(:,burnin_ite+autocorre_time*s);
  
end
fprintf('label switching...\n')
   latent_vector=labelswitch(z_pos);  % label switching
   %latent_vector=z_pos;
   % Estimate z (the most frequent value)
   latent_pos=latent_vector';  
   [Au,ia,ic] = unique(latent_pos,'rows','stable');
   Counts = accumarray(ic, 1);
   Out = [Counts Au];
   output=sortrows(Out,1,'descend'); % sort in descend
   z_esti=output(1,2:end)'; % the first row is the most frequent
   % Estimate K (the most frequent value)
   K=mode(K_pos);
end

 
%--------------------------------------------------------------------------
function [z_chain,prob_chain,accep_r,K_chain]=MCMC_allocation_sampler(x,K,priorpra)

% Metroplis-Hastings with Gibbs move, M3 move and AE move
%
% Input: 
%   x: the observation, a correlation matrix of the network
%   K: Initialize the number of communities
%   priorpra: a struct containing the prior hyper parameters
%
% Output: z_chain: a Markov chain, N x Itera matrix
%         prob_chain: posterior probability of the state
%         accep_r: acceptance ratio of each step
%         K_chain: number of communities of each state
%
% Version 1.0
% 15-March-2020
% Copyright (c) 2020, Lingbin Bian
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Markov chain Monte Carlo

a=1;

N=length(x);
latent_ini=latent_initial(N,K);   % initialization

Itera=priorpra.Itera;        % number of iteration
z_chain=zeros(N,Itera);
prob_chain=zeros(1,Itera);
accep_r=zeros(1,Itera);
K_chain=zeros(1,Itera);

z_chain(:,1)=latent_ini;
prob_chain(1,1)=log(1/(factorial(K))*exp(-1))...
                +pos_z_K(latent_ini,K,N)+gauss_x_z(latent_ini,x,K,N);
K_chain(1,1)=max(latent_ini);


for ite=2:Itera
   
    fprintf('MCMC iteration, total steps: 2000, step: %d\n',ite);
    p_AE=0.5;
    ind=binornd(1,p_AE);
    %ind=randi(2);
    if ind==0
        p_gibbs=0.5;
        ind_gibbs=binornd(1,p_gibbs);
        if ind_gibbs==0 && K_chain(1,ite-1)>1 
          [z_chain(:,ite),prob_chain(1,ite),accep_r(1,ite)]...
              =move_M_3(z_chain(:,ite-1),x,K_chain(1,ite-1),N);
           K_chain(1,ite)=max(z_chain(:,ite));           
        else
            [z_chain(:,ite),prob_chain(1,ite),accep_r(1,ite)]...
              =move_gibbs(z_chain(:,ite-1),x,K_chain(1,ite-1),N);
           K_chain(1,ite)=max(z_chain(:,ite)); 
        end            
    elseif ind==1
       [z_chain(:,ite),prob_chain(1,ite),accep_r(1,ite),K_chain(1,ite)]...
           =move_AE(z_chain(:,ite-1),x,N);
    end        
end

% -------------------------------------------------------------------------
    % nested functions for function [z_chain,prob_chain,accep_r,K_chain]=MCMC_allocation_sampler(x,K,priorpra)
    %----------------------------------------------------------------------
    % Initialize latent label vector     
     function [z] = latent_initial(N,K)
        % This function initializes the latent label vector by Uniform distribution.
        leng=0;
        while leng~=K
        z=randi([1,K],N,1);
        leng=length(unique(z));
        end
     end


     % Gibbs sampling

     function [ z, prob_chain,accep_r ] = move_gibbs(z,x,K,N )  % Computational Complexity: O((K^2+N^2)K)
            % p(z_i|z_-i,x)  for i=1:N
            % a=log(p_1),b=log(p_2),log_sum=log(p_1+p_2)    
            prob=zeros(N,K);
            prob_n=zeros(N,K);
            log_sum=zeros(N,1);
            i=randi(N);
            
            for k=1:K
                z(i)=k;            
                prob(i,k)=log(1/(factorial(K))*exp(-1))+pos_z_K(z,K,N)+gauss_x_z(z,x,K,N);            
                %               prob(i,k)=log(1/(factorial(K))*exp(-1))+Pos_z_K(z,K,N)+Gauss_x_z(z,x,K,N) % log(p(x|z]))
                %               prob_sum(i)=prob_sum(i)+prob(i,k);
            end
            %           log_sum(i)=prob(i,1);
            %           for k=2:K
            %               log_sum(i)=ln_sum(log_sum(i),prob(i,k));
            %           end
            log_sum(i)=logsumexp(prob(i,:)');
            for k=1:K
                prob_n(i,k)=prob(i,k)-log_sum(i);% normalized posterior probability            
                prob_n(i,k)=exp(prob_n(i,k));          
            end
            z_vec=mnrnd(1,single(prob_n(i,:)),1); % latent label of vector form
            z(i)=find(z_vec==1); % transfer from vector form to scaler eg: [0 0 1 0] to 3   
            prob_chain=log(1/(factorial(K))*exp(-1))+pos_z_K(z,K,N)+gauss_x_z(z,x,K,N);
            accep_r=1;
     end
        
    
    %----------------------------------------------------------------------
    
    %ln_sum=@(a,b) a+log(1+exp(b-a));
    
        function [s] = logsumexp(a, dim)  % O(1)
        % Returns log(sum(exp(a),dim)) while avoiding numerical underflow.
        % Default is dim = 1 (columns).
        % logsumexp(a, 2) will sum across rows instead of columns.
        % Unlike matlab's "sum", it will not switch the summing direction
        % if you provide a row vector.
         
        % Written by Tom Minka
        % (c) Microsoft Corporation. All rights reserved.
    
            if nargin < 2
              dim = 1;
            end
    
            % subtract the largest in each column
            [y, i] = max(a,[],dim);
            dims = ones(1,ndims(a));
            dims(dim) = size(a,dim);
            a = a - repmat(y, dims);
            s = y + log(sum(exp(a),dim));
            i = find(~isfinite(y));
            if ~isempty(i)
              s(i) = y(i);
            end
        end
    
    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------
    % MCMC with M3
    function [z_star,prob_chain,accep_r]=move_M_3(z,x,K,N)   % K^2+N^2+L^2
    % This function is the Absorption-Ejection move
    % Input: z: current state
    %        x: observation, a correlation matrix
    %        K: maximum value of community number
    %        N: number of nodes
    % Output: z_star: the updated state
    %         prob_chain: posterior probability
    %         accep_r: acceptance ratio
        
        [z_update,propmove_ratio]=proposal_move(z,x,K,N);   % M3 move 
        r=pos(z_update,x,K,N)-pos(z,x,K,N)+log(propmove_ratio);
        r=exp(r);
        accep_r=min(1,r);
        u=rand();    
        if u<accep_r
            z_star=z_update;
        else
            z_star=z;
        end
        prob_chain=log(1/(factorial(K))*exp(-1))+pos_z_K(z,K,N)...
                   +gauss_x_z(z,x,K,N);
    end
        
    % ---------------------------------------------------------------------
        function [z_update,propmove_ratio]=proposal_move(z,x,K,N) % Computational cost: N+L^2
        % The M3 move 
        % randomly select two clusters
        
        p=randsample(1:K,2);
        k_1=p(1);  % k1
        k_2=p(2);  % k2
        % number of nodes in k1 and k2
        m_k=zeros(K,1); %       
        for i=1:N
            m_k(z(i))=m_k(z(i))+1;
        end
        Lk_1=m_k(k_1);   % m_k_1
        Lk_2=m_k(k_2);   % m_k_2
    
        % whole label vector [index,labels]
        z_label=zeros(N,2);
        z_label(:,1)=1:N;
        z_label(:,2)=z;   
        % list in A={i: z_i = k1 or k2}
        z_list=zeros(Lk_1+Lk_2,2);
        % list not in A
        z_tilde=zeros(N-Lk_1-Lk_2,2);
        
        n_tilde=1;
        n_list=1;
        for i=1:N
            if z_label(i,2)==k_1 || z_label(i,2)==k_2
                z_list(n_list,:)=z_label(i,:);
                n_list=n_list+1;       
            else
                z_tilde(n_tilde,:)=z_label(i,:);
                n_tilde=n_tilde+1;
            end
        end
      
        L_list=length(z_list);    
        propmove_ratio=1;
        while L_list~=0    % L(N-L+L^2)
            [z_tilde,z_list,L_list,propmove_ratio]...
                =z_tilde_update(x,z_tilde,z_list,propmove_ratio,k_1,k_2); 
        end
        z_tilde=sortrows(z_tilde,1);
        z_update=z_tilde(:,2);
        end
    
    %----------------------------------------------------------------------
        function [z_tilde,z_list,L_list,propmove_ratio]...
                =z_tilde_update(x,z_tilde,z_list,propmove_ratio,k_1,k_2)  % N-L+L^2
        % This function updates z~ and the list of z in A={i: z_i = k1 or k2}
           
            L_list=length(z_list(:,1));  % length of list
            l=randi(L_list);    % randomly select a seed in the list
            z_seed=z_list(l,:);
            L_tilde=length(z_tilde(:,1));  % length of z~
    
            m_k_1=0;
            m_k_2=0;
            for i=1:L_tilde
                if z_tilde(i,2)==k_1
                    m_k_1=m_k_1+1;   % m_k1 tilde
                elseif z_tilde(i,2)==k_2
                    m_k_2=m_k_2+1;   % m_k2 tilde
                end
            end
            z_tilde_k_1=zeros(m_k_1,2);
            z_tilde_k_2=zeros(m_k_2,2);
    
            l_1=1;
            l_2=1;
    
            for i=1:L_tilde
                if z_tilde(i,2)==k_1
                    z_tilde_k_1(l_1,:)=z_tilde(i,:);
                    l_1=l_1+1;
    
                elseif z_tilde(i,2)==k_2
                    z_tilde_k_2(l_2,:)=z_tilde(i,:);
                    l_2=l_2+1;
                end
            end
    
            z_tilde_k_1_ws=[z_tilde_k_1(:,:);z_seed(1,1),k_1];
            z_tilde_k_2_ws=[z_tilde_k_2(:,:);z_seed(1,1),k_2];
            %z_tilde_k_1_ws=sortrows(z_tilde_k_1_ws);
            %z_tilde_k_2_ws=sortrows(z_tilde_k_2_ws);
    
            r_1=prob_ratio_1(m_k_1,m_k_2);
    
            r_2=exp(prob_ratio_2(x,z_tilde_k_1,z_tilde_k_2,...
                z_tilde_k_1_ws,z_tilde_k_2_ws,m_k_1,m_k_2));
    
            r=r_1*r_2;
            prob_k_1=r/(1+r);
            prob_k_2=1-prob_k_1;
    
            % reverse assignment probability
            if z_seed(1,2)==k_1
                prop_reverse=prob_k_1;
            else
                prop_reverse=prob_k_2;
            end
    
            % update seed and forward assignment probability
            Ber_num=binornd(1,prob_k_1);
            if Ber_num==1
                z_seed(1,2)=k_1;
                prop_forward=prob_k_1;
            else
                z_seed(1,2)=k_2;
                prop_forward=prob_k_2;
            end
            propmove_ratio=propmove_ratio*prop_reverse/prop_forward;
    
            z_tilde=[z_tilde(:,:);z_seed(1,:)];
            %z_tilde=sortrows(z_tilde,1);
    
            z_list(l,:)=[];
            L_list=L_list-1;
        end
    %----------------------------------------------------------------------
    % r_1=p(z*_i=k_1,z~|K)/p(z*_i=k_2,z~|K)
        function [r_1]=prob_ratio_1(num_k_1,num_k_2)
            r_1=(1+num_k_1)/(1+num_k_2);        
        end
    %----------------------------------------------------------------------
    % r_2=p(x*_i,x~|K,z*_i=k1,z~)/p(x*_i,x~|K,z*_i=k2,z~)
        function [r_2]=prob_ratio_2(x,z_tilde_k_1,z_tilde_k_2,...   % 
                       z_tilde_k_1_ws,z_tilde_k_2_ws,m_k_1,m_k_2)
            if m_k_1>0&&m_k_2>0         
                r_2=p_b(x(z_tilde_k_1_ws(:,1),z_tilde_k_1_ws(:,1)))...
                   +p_b(x(z_tilde_k_2(:,1),z_tilde_k_2(:,1)))...
                   -p_b(x(z_tilde_k_1(:,1),z_tilde_k_1(:,1)))...
                   -p_b(x(z_tilde_k_2_ws(:,1),z_tilde_k_2_ws(:,1)));
            elseif m_k_1>0&&m_k_2==0
                r_2=p_b(x(z_tilde_k_1_ws(:,1),z_tilde_k_1_ws(:,1)))...
                    -p_b(x(z_tilde_k_1(:,1),z_tilde_k_1(:,1)))...
                    -p_b(x(z_tilde_k_2_ws(:,1),z_tilde_k_2_ws(:,1)));
            elseif m_k_1==0&&m_k_2>0  
                r_2=p_b(x(z_tilde_k_1_ws(:,1),z_tilde_k_1_ws(:,1)))...
                    +p_b(x(z_tilde_k_2(:,1),z_tilde_k_2(:,1)))...
                    -p_b(x(z_tilde_k_2_ws(:,1),z_tilde_k_2_ws(:,1)));
            elseif m_k_1==0&&m_k_2==0
                r_2=p_b(x(z_tilde_k_1_ws(:,1),z_tilde_k_1_ws(:,1)))...
                    -p_b(x(z_tilde_k_2_ws(:,1),z_tilde_k_2_ws(:,1)));        
            end
        end
    %----------------------------------------------------------------------
        function [p_block]=p_b(x_block) 
        % block log likelihood
             nu=priorpra.nu;    % sigma_kl^2 ~ IG(nu/2,rho/2)
             rho=priorpra.rho;    
             xi=priorpra.xi;   % mu_kl ~ N(xi,kappa^2*sigma^2)
             kappa_sq=priorpra.kappa_sq;
             N_b=length(x_block);
                   
             W_b=N_b^2;
             w_sum=0;
             w_sumsq=0;
             for i=1:N_b
                 for j=1:N_b             
                     w_sum=w_sum+x_block(i,j);
                     w_sumsq=w_sumsq+x_block(i,j).^2;                
                 end
             end
             Term_A=(0.5*nu)*log(rho);        
             Term_B=gammaln(0.5*(W_b+nu));
             Term_C=(0.5*W_b)*log(pi);
             Term_D=gammaln(0.5*nu);
             Term_E=log((W_b*kappa_sq+1)^0.5);
             Term_F=(-0.5*(W_b+nu))...
                 *log((w_sumsq-(kappa_sq*((w_sum+xi/kappa_sq)^2))...
                 /(W_b*kappa_sq+1)+xi^2/kappa_sq+rho));
             p_block=Term_A+Term_B-Term_C-Term_D-Term_E+Term_F;
        end
    %----------------------------------------------------------------------
    % AE move -------------------------------------------------------------
    function [z_star,prob_chain,accep_r,K_star]=move_AE(z,x,N)
    % This function is the absorption-ejection move
    % Input: z: current state
    %        x: observation, a correlation matrix
    %        N: number of nodes
    % Output: z_star: the updated state
    %         prob_chain: posterior probability
    %         accep_r: acceptance ratio
    
        K_max=20;
        if max(z)==K_max
            p_ek=0;
        elseif max(z)==1
            p_ek=1;
        else
            p_ek=0.5;
        end
        AE_ind=binornd(1,p_ek);
        if AE_ind==1
            [z_update,propmove_ratio]=proposal_eject(z,N);
        elseif AE_ind==0
            [z_update,propmove_ratio]=proposal_absorb(z,N);
        end        
        r=pos(z_update,x,max(z_update),N)-pos(z,x,max(z),N)+log(propmove_ratio);
        r=exp(r);
        accep_r=min(1,r);
        u=rand();
        if u<accep_r
            z_star=z_update;
        else
            z_star=z;
        end
        prob_chain=log(1/(factorial(max(z_star)))*exp(-1))+pos_z_K(z_star,max(z_star),N)...
            +gauss_x_z(z_star,x,max(z_star),N);
        K_star=max(z_star);
    end
            
    %----------------------------------------------------------------------
    % ejection
            function [z_update,propmove_ratio]=proposal_eject(z,N)             
                K_z=max(z);
                K_r=randsample(1:K_z,1);           
                k_1=K_r;
                k_2=K_z+1;  
                P_E=betarnd(double(a),double(a));
                
                for i=1:N
                    if z(i)==k_1
                        Ber_num=binornd(1,P_E);
                        if Ber_num==1
                            z(i)=k_2;
                        else
                            z(i)=k_1;
                        end
                    end
                end            
                z_update=z;            
                m_k=zeros(K_z+1,1); %
                for i=1:N
                    m_k(z(i))=m_k(z(i))+1;
                end
                n_k_1=m_k(k_1);
                n_k_2=m_k(k_2);
                n=n_k_1+n_k_2;
                propmove_ratio=prop_ratio(n_k_1,n_k_2,n);            
            end
    %----------------------------------------------------------------------  
    % absorption
            function [z_update,propmove_ratio]=proposal_absorb(z,N)           
                K_z=max(z);
                K_r=randsample(1:(K_z-1),1);
                k_1=K_r;
                k_2=K_z;            
                m_k=zeros(K_z,1); %
                for i=1:N
                    m_k(z(i))=m_k(z(i))+1;
                end
                n_k_1=m_k(k_1);
                n_k_2=m_k(k_2);
                n=n_k_1+n_k_2;
                for i=1:N
                    if z(i)==k_2
                        z(i)=k_1;
                    end
                end
                z_update=z;            
                propmove_ratio=1/prop_ratio(n_k_1,n_k_2,n);            
            end   
    %----------------------------------------------------------------------
    % proposal ratio
            function [propmove_ratio]=prop_ratio(n_k_1,n_k_2,n)            
                propmove_ratio=gamma(a)*gamma(a)/(gamma(2*a))*(gamma(2*a+n)/(gamma(a+n_k_1)*gamma(a+n_k_2)));            
            end
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % Posterior probability -----------------------------------------------
    
        function [posteri]=pos(z,x,K,N)
            posteri=log(1/(factorial(K))*exp(-1))...
                    +pos_z_K(z,K,N)+gauss_x_z(z,x,K,N);
        end
    %----------------------------------------------------------------------
    % log(p(x|z))
        function [ P_x_z] = gauss_x_z(z,x,K,N)   % O(K^2+N^2)
            W=zeros(K,K);     % number of elements in the block
            w_sum=zeros(K,K);    % weights in the block
            w_sumsq=zeros(K,K);   % weight^2 in the block
            m_k=zeros(K,1);   % number of nodes in the cluster     
    
            nu=priorpra.nu;    % sigma_kl^2 ~ IG(nu/2,rho/2)
            rho=priorpra.rho;    
            xi=priorpra.xi;   % mu_kl ~ N(xi,kappa^2*sigma^2)
            kappa_sq=priorpra.kappa_sq;
    
            for i=1:N
                m_k(z(i))=m_k(z(i))+1;
            end      
            for k=1:K
                for l=1:K               
                   W(k,l)=m_k(k)*m_k(l);    % w_kl=m_k*m_l        
                end
            end      
            for i=1:N
                for j=1:N     
                    w_sum(z(i),z(j))=w_sum(z(i),z(j))+x(i,j);
                    w_sumsq(z(i),z(j))=w_sumsq(z(i),z(j))+x(i,j).^2;              
                end
            end 
    
            P_x_z=zeros(K,K);
            for k=1:K
                for l=1:K   
                    Term_A=(0.5*nu)*log(rho);        
                    Term_B=gammaln(0.5*(W(k,l)+nu));
                    Term_C=(0.5*W(k,l))*log(pi);
                    Term_D=gammaln(0.5*nu);
                    Term_E=log((W(k,l)*kappa_sq+1)^0.5);
                    Term_F=(-0.5*(W(k,l)+nu))...
                        *log((w_sumsq(k,l)...
                        -(kappa_sq*((w_sum(k,l)+xi/kappa_sq)^2))...
                        /(W(k,l)*kappa_sq+1)+xi^2/kappa_sq+rho));
                    P_x_z(k,l)=Term_A+Term_B-Term_C-Term_D-Term_E+Term_F;
    
        % P_x_z(k,l)=log(rho^(0.5*nu))+gammaln(0.5*(W(k,l)+nu))...
        % -log(pi^(0.5*W(k,l)))-gammaln(0.5*nu)-log((W(k,l)*kappa_sq+1)^0.5)...
        % +log((w_sumsq(k,l)-(kappa_sq*((w_sum(k,l)+xi/kappa_sq)^2))...
        % /(W(k,l)*kappa_sq+1)+xi^2/kappa_sq+rho)^(-0.5*(W(k,l)+nu)));  
                end
            end
    
            P_x_z=sum(sum(P_x_z)); % log(p(x|z))
        end
    %----------------------------------------------------------------------
    % log(p(z|K))
        function [ P_z_K ] = pos_z_K(z,K,N)  % O(N+K)
            m_k=zeros(K,1); % number of nodes in the cluster
          
            for i=1:N
                m_k(z(i))=m_k(z(i))+1;
            end
    %         for k=1:K
    %             m_k(k)=length(find(z==k));
    %         end
            Gam_1M=zeros(K,1); 
            for i=1:K
                Gam_1M(i)=gammaln(m_k(i)+1)-gammaln(1);
            end
            term=sum(Gam_1M);
            P_z_K=gammaln(K)-gammaln(K+N)+term;
        end
    
    % End of nested functions----------------------------------------------
end
 


%--------------------------------------------------------------------------
% label switching

function [label_v]=labelswitch(label_v)
% This function solves the label switching phenomenon.
% label_v: N*M matrix, each column is a latent label vector
%          N is the number of nodes, 
%          M is the number of latent lable vectors
%
% Version 1.0
% 7-Aug-2020
% Copyright (c) 2020, Lingbin Bian

% label_v=sort_vect(label_v);
[n_row,n_column]=size(label_v);
cost_mat=cell(1,n_column);
sigma=cell(1,n_column);
for j=2:n_column    
   cost_mat{1,j}=zeros(max(label_v(:,j)),max(label_v(:,j)));
   for k=1:max(label_v(:,j))
       for l=1:max(label_v(:,j))
           for t=1:j-1
               for i=1:n_row
                   if label_v(i,j-t)~=k && label_v(i,j)==l
                       cost_mat{1,j}(k,l)=cost_mat{1,j}(k,l)+1;
                   end
               end
           end
       end
   end

   sigma{1,j}=munkres(cost_mat{1,j});
   L_sig=length(sigma{1,j});
   for i=1:n_row
       inde=0;
       for k_sigma=1:L_sig
       if label_v(i,j)==sigma{1,j}(k_sigma)
           label_v(i,j)=k_sigma;
           inde=1;
       end
       if inde==1        
           break;           
       end
       end           
   end      
end
             
end


function [assignment,cost] = munkres(costMat)
        % MUNKRES   Munkres (Hungarian) Algorithm for Linear Assignment Problem. 
        %
        % [ASSIGN,COST] = munkres(COSTMAT) returns the optimal column indices,
        % ASSIGN assigned to each row and the minimum COST based on the assignment
        % problem represented by the COSTMAT, where the (i,j)th element represents the cost to assign the jth
        % job to the ith worker.
        %
        % Partial assignment: This code can identify a partial assignment is a full
        % assignment is not feasible. For a partial assignment, there are some
        % zero elements in the returning assignment vector, which indicate
        % un-assigned tasks. The cost returned only contains the cost of partially
        % assigned tasks.
        
        % This is vectorized implementation of the algorithm. It is the fastest
        % among all Matlab implementations of the algorithm.
        
        % Examples
        % Example 1: a 5 x 5 example
        %{
        [assignment,cost] = munkres(magic(5));
        disp(assignment); % 3 2 1 5 4
        disp(cost); %15
        %}
        % Example 2: 400 x 400 random data
        %{
        n=400;
        A=rand(n);
        tic
        [a,b]=munkres(A);
        toc                 % about 2 seconds 
        %}
        % Example 3: rectangular assignment with inf costs
        %{
        A=rand(10,7);
        A(A>0.7)=Inf;
        [a,b]=munkres(A);
        %}
        % Example 4: an example of partial assignment
        %{
        A = [1 3 Inf; Inf Inf 5; Inf Inf 0.5]; 
        [a,b]=munkres(A)
        %}
        % a = [1 0 3]
        % b = 1.5
        % Reference:
        % "Munkres' Assignment Algorithm, Modified for Rectangular Matrices", 
        % http://csclab.murraystate.edu/bob.pilgrim/445/munkres.html
        
        % version 2.3 by Yi Cao at Cranfield University on 11th September 2011
        
        assignment = zeros(1,size(costMat,1));
        cost = 0;
        
        validMat = costMat == costMat & costMat < Inf;
        bigM = 10^(ceil(log10(sum(costMat(validMat))))+1);
        costMat(~validMat) = bigM;
        
        % costMat(costMat~=costMat)=Inf;
        % validMat = costMat<Inf;
        validCol = any(validMat,1);
        validRow = any(validMat,2);
        
        nRows = sum(validRow);
        nCols = sum(validCol);
        n = max(nRows,nCols);
        if ~n
            return
        end
        
        maxv=10*max(costMat(validMat));
        
        dMat = zeros(n) + maxv;
        dMat(1:nRows,1:nCols) = costMat(validRow,validCol);
        
        %*************************************************
        % Munkres' Assignment Algorithm starts here
        %*************************************************
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   STEP 1: Subtract the row minimum from each row.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        minR = min(dMat,[],2);
        minC = min(bsxfun(@minus, dMat, minR));
        
        %**************************************************************************  
        %   STEP 2: Find a zero of dMat. If there are no starred zeros in its
        %           column or row start the zero. Repeat for each zero
        %**************************************************************************
        zP = dMat == bsxfun(@plus, minC, minR);
        
        starZ = zeros(n,1);
        while any(zP(:))
            [r,c]=find(zP,1);
            starZ(r)=c;
            zP(r,:)=false;
            zP(:,c)=false;
        end
        
        while 1
        %**************************************************************************
        %   STEP 3: Cover each column with a starred zero. If all the columns are
        %           covered then the matching is maximum
        %**************************************************************************
            if all(starZ>0)
                break
            end
            coverColumn = false(1,n);
            coverColumn(starZ(starZ>0))=true;
            coverRow = false(n,1);
            primeZ = zeros(n,1);
            [rIdx, cIdx] = find(dMat(~coverRow,~coverColumn)==bsxfun(@plus,minR(~coverRow),minC(~coverColumn)));
            while 1
                %**************************************************************************
                %   STEP 4: Find a noncovered zero and prime it.  If there is no starred
                %           zero in the row containing this primed zero, Go to Step 5.  
                %           Otherwise, cover this row and uncover the column containing 
                %           the starred zero. Continue in this manner until there are no 
                %           uncovered zeros left. Save the smallest uncovered value and 
                %           Go to Step 6.
                %**************************************************************************
                cR = find(~coverRow);
                cC = find(~coverColumn);
                rIdx = cR(rIdx);
                cIdx = cC(cIdx);
                Step = 6;
                while ~isempty(cIdx)
                    uZr = rIdx(1);
                    uZc = cIdx(1);
                    primeZ(uZr) = uZc;
                    stz = starZ(uZr);
                    if ~stz
                        Step = 5;
                        break;
                    end
                    coverRow(uZr) = true;
                    coverColumn(stz) = false;
                    z = rIdx==uZr;
                    rIdx(z) = [];
                    cIdx(z) = [];
                    cR = find(~coverRow);
                    z = dMat(~coverRow,stz) == minR(~coverRow) + minC(stz);
                    rIdx = [rIdx(:);cR(z)];
                    cIdx = [cIdx(:);stz(ones(sum(z),1))];
                end
                if Step == 6
                    % *************************************************************************
                    % STEP 6: Add the minimum uncovered value to every element of each covered
                    %         row, and subtract it from every element of each uncovered column.
                    %         Return to Step 4 without altering any stars, primes, or covered lines.
                    %**************************************************************************
                    [minval,rIdx,cIdx]=outerplus(dMat(~coverRow,~coverColumn),minR(~coverRow),minC(~coverColumn));            
                    minC(~coverColumn) = minC(~coverColumn) + minval;
                    minR(coverRow) = minR(coverRow) - minval;
                else
                    break
                end
            end
            %**************************************************************************
            % STEP 5:
            %  Construct a series of alternating primed and starred zeros as
            %  follows:
            %  Let Z0 represent the uncovered primed zero found in Step 4.
            %  Let Z1 denote the starred zero in the column of Z0 (if any).
            %  Let Z2 denote the primed zero in the row of Z1 (there will always
            %  be one).  Continue until the series terminates at a primed zero
            %  that has no starred zero in its column.  Unstar each starred
            %  zero of the series, star each primed zero of the series, erase
            %  all primes and uncover every line in the matrix.  Return to Step 3.
            %**************************************************************************
            rowZ1 = find(starZ==uZc);
            starZ(uZr)=uZc;
            while rowZ1>0
                starZ(rowZ1)=0;
                uZc = primeZ(rowZ1);
                uZr = rowZ1;
                rowZ1 = find(starZ==uZc);
                starZ(uZr)=uZc;
            end
        end
        
        % Cost of assignment
        rowIdx = find(validRow);
        colIdx = find(validCol);
        starZ = starZ(1:nRows);
        vIdx = starZ <= nCols;
        assignment(rowIdx(vIdx)) = colIdx(starZ(vIdx));
        pass = assignment(assignment>0);
        pass(~diag(validMat(assignment>0,pass))) = 0;
        assignment(assignment>0) = pass;
        cost = trace(costMat(assignment>0,assignment(assignment>0)));
        
        function [minval,rIdx,cIdx]=outerplus(M,x,y)
        ny=size(M,2);
        minval=inf;
        for c=1:ny
            M(:,c)=M(:,c)-(x+y(c));
            minval = min(minval,min(M(:,c)));
        end
        [rIdx,cIdx]=find(M==minval);
        end
        
end
        % end of 2nd layer nested functions