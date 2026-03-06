% Community detection using modularity (synthetic data)
% individual-level analysis
%
% Version 1.0
% 23-Feb-2026
% Copyright (c) 2026, Lingbin Bian
% -------------------------------------------------------------------------
clear
clc
close all
% -------------------------------------------------------------------------
% Load synthetic data
datatype=0;
subjects=load('synthetic_id.txt'); 
session_n='';

% signal to noise ratio SNR=10log10(sigma_signal^2/sigma_noise^2);
% sigma_signal=1, sigma_noise=n_s
n_s=0.3162;  % 10dB
% n_s=0.5623;  % 5dB
% n_s=1;      
% n_s=1.7783;  % -5dB
% n_s=3.1623;  % -10dB
% n_s=5.6234;  % -15dB
% n_s=10;      % -20dB
% n_s=17.7828; % -25dB

atlas=1;
DIIV=[10,20,30]; 

opt_para_modularity=zeros(3,3);
for state=1:3
    for n_diiv=1:length(DIIV)
        load(['Results/synthetic_Bayes_opt_M','/','DIIV',num2str(DIIV(n_diiv)),'/n',num2str(n_s),'/','opt_para_',num2str(state)]);
       
        opt_para_modularity(state,n_diiv)=bestGamma;
    end
end


hrf_ind=0;

state_t=[50,150,250]; 

K_state=[8,9,10];
W=10;

N_subj=100;

N_state=length(state_t);

group_adj=cell(N_subj,length(DIIV)); % group of adjacency matrix subject x DIIV

Q=cell(N_subj,length(DIIV));
K_g=cell(N_subj,length(DIIV)); % group of estimated K
z_g=cell(N_subj,length(DIIV)); % group of estimated z
true_latent_labels=cell(1,length(DIIV));

nmi_score=cell(N_state,length(DIIV));

for ds=1:N_state
    for i=1:length(DIIV)
        nmi_score{ds,i}=zeros(N_subj,1);
    end
end

for s=1:N_subj
    fprintf('Adjacency matrix of subject: %d\n',s)
    subid=num2str(subjects(s));
    for i=1:length(DIIV)
    [group_adj{s,i},true_latent,true_latent_labels{1,i}]=local_adj(datatype,atlas,subid,session_n,n_s,state_t,K_state,W,DIIV(i),hrf_ind);
    end
end

% -------------------------------------------------------------------------
% Adjacency matrix

M=10; % margin size

N=100; 

for s=1:N_subj
    fprintf('subject: %d\n',s)
    for i=1:length(DIIV)
        z_ds=zeros(N,N_state);
        Q_ds=zeros(1,N_state);
        K_ds=zeros(1,N_state);
        for ds=1:N_state
          [z_ds(:,ds),Q_ds(1,ds)] = modularity_und(group_adj{s,i}{1,ds},opt_para_modularity(ds,i));   
          K_ds(1,ds)=max(z_ds(:,ds));
        end
        z_g{s,i}=z_ds;
        Q{s,i}=Q_ds;
        K_g{s,i}=K_ds;
    end
end

for ds=1:N_state
    for i=1:length(DIIV)        
        for s=1:N_subj
            nmi_score{ds,i}(s,:)=nmi(z_g{s,i}(:,ds),true_latent_labels{1,i}{s,1}(:,ds));
        end        
    end
end

% Save the results of individual-level modelling

data_path = fileparts(mfilename('fullpath'));

group_path=fullfile(data_path,'Results/synthetic_Bayes_opt_M/infer_opt');
save(group_path,'opt_para_modularity','z_g','true_latent_labels','nmi_score');
   












