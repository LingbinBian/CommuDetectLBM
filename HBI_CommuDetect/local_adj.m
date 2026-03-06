function [adj,true_latent,true_latent_sub,K_seg]=local_adj(datatype,atlas,subj,session_n,n_s,state_t,K_state,W,vari,hrf_ind)
% This function calculates the local adjacency matrix of single subject.
% 
% Input: datatype: 1 real, 2 synthetic 
%        subj: subject ID
%        session_n: 1: LR, 2 RL (only for real data)
%        localmin_t: a vector of time points with respect to the local minima
%        of the time series
%        K_min: a vector containing number of communities of the states
%        W: half of window size
% Output: adj: local adjacency matrix
%         true_latent: true community memberships
%         K_seg: true number of communities
%
% Version 1.0
% 11-Jun-2020
% Copyright (c) 2020, Lingbin Bian
% -------------------------------------------------------------------------
% load data

if datatype==1
    Timeseries=load_realts(atlas,subj,session_n,K_state);
    true_latent='none';
    true_latent_sub='none';
    K_seg=K_state;
elseif datatype==0
    [Timeseries,latent_seg,latent_seg_sub,K_seg]=load_simults(subj,n_s,vari,hrf_ind);
    true_latent=latent_seg;
    true_latent_sub=latent_seg_sub;
end

L_state_t=length(state_t);   % number of states

adj=cell(1,L_state_t);  % each cell is an adjacency matrix

for j=1:L_state_t
    adj{j}=adjac_generator(state_t(j));
end


% Nested functions---------------------------------------------------------
% This function generates the adjacency matrix within the window at time t
    function [ Adj] = adjac_generator(t)
        signal_inW=Timeseries.signal(:,t-W:t+W-1)';
        covari=cov(signal_inW);
        Adj=corrcov(covari);
%        Adj=corrcoef(signal_inW);
%        Adj=threshold_absolute(corre,0.2);
%         for i=1:N
%             for j=1:N
%                 if Adj(i,j)>=0.2
%                     Adj(i,j)=1;
%                 else
%                     Adj(i,j)=0;
%                 end
%             end
%         end        
    end
end

