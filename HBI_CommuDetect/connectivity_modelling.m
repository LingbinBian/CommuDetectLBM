function [esti_connectmean,esti_connectvariance]=connectivity_modelling(group_adj,N,localmin_t,S)
% This function estimates block mean and variance.
% 
% Input: adj_state: adjacency cell for all of the states
%        localmin_t: a vector of time points with respect to the states
%        K_min: a vector containing number of communities of states
%        S: replication number
% Output: esti_blockmean: estimation of connectivity mean
%         esti_variance: estimation of connectivity variance
%
% Version 1.0
% 11-April-2021
% Copyright (c) 2021, Lingbin Bian
% -------------------------------------------------------------------------
 
S_subj=100;
L_localmin=length(localmin_t);   % number of states

adj_state=cell(1,L_localmin);

% Convert group_adj to adj_statef
for s_st=1:L_localmin
   adj_state{1,s_st}=cell(N,N);
   for i=1:N
       for j=1:N        
            X_ij=zeros(S_subj,1);
            for s_j=1:S_subj
                 X_ij(s_j,1)=group_adj{s_j,1}{1,s_st}(i,j);
            end
            adj_state{1,s_st}{i,j}=X_ij;                 
       end
   end
end

L_localmin=length(localmin_t);   % number of states

connect_parameter=cell(1,L_localmin);  % each cell is a struct containing connectivity parameters
cell_connectmean=cell(S,L_localmin);  % a cell containing replications of block mean
cell_connectvariance=cell(S,L_localmin);  % a cell containing replications of block variance
esti_connectmean=cell(1,L_localmin);  % estimated connectivity mean at each state
esti_connectvariance=cell(1,L_localmin);  % estimated connectivity variance at each state

for j=1:L_localmin
    connect_parameter{1,j}=element_modelling(adj_state{j},N,S);
end

for j=1:L_localmin
   for s=1:S
       para=connect_parameter{1,j}(s);
       cell_connectmean{s,j}=para.Mean;
       cell_connectvariance{s,j}=para.Vari;
   end
end

Mean_ij=zeros(S,1); % store temporal single mean element of replications 
for t=1:L_localmin
    for i=1:N
        for j=1:N
            for s=1:S
                Mean_ij(s,1)=cell_connectmean{s,t}(i,j);
            end
            esti_connectmean{1,t}(i,j)=mean(Mean_ij);
        end
    end
end

Variance_ij=zeros(S,1); % store temporal single variance element of replications
for t=1:L_localmin
    for i=1:N
        for j=1:N
            for s=1:S
                Variance_ij(s,1)=cell_connectvariance{s,t}(i,j);
            end
            esti_connectvariance{1,t}(i,j)=mean(Variance_ij);
        end
    end
end

end


