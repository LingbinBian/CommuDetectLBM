% Community detection based on multilayer modularity
% synthetic data analysis
%
% Version 1.0
% 1-March-2026
% Copyright (c) 2026, Lingbin Bian
% -------------------------------------------------------------------------
clear
clc
close all
% -------------------------------------------------------------------------
% Load data
datatype=0;
subjects=load('synthetic_id.txt'); 
session_n='';
atlas=1;
% signal to noise ratio SNR=10log10(sigma_signal^2/sigma_noise^2);
% sigma_signal=1, sigma_noise=n_s
n_s=0.3162;  % 10dB
% n_s=0.5623;  % 5d
% n_s=1;      
% n_s=1.7783;  % -5dB
% n_s=3.1623;  % -10dB
% n_s=5.6234;  % -15dB
% n_s=10;      % -20dB
% n_s=17.7828; % -25dB


DIIV=[10,20,30];

opt_para_multi_modularity=cell(3,3);
for state=1:3
    for n_diiv=1:length(DIIV)
        opt_para_multi_modularity{state,n_diiv}=zeros(2,1);
    end
end

% load the optimal parameters
for state=1:3
    for n_diiv=1:length(DIIV)
    load(['Results/synthetic_Bayes_opt_MM','/','DIIV',num2str(DIIV(n_diiv)),'/n',num2str(n_s),'/','opt_para_',num2str(state)]);
    opt_para_multi_modularity{state,n_diiv}(1,1)=bestGamma; %K_max
    opt_para_multi_modularity{state,n_diiv}(2,1)=bestOmega; %nu
    end
end


para=opt_para_multi_modularity;
% coupling parameter

vari=[10,20,30]; % true community structure
hrf_ind=0;

state_t=[50,150,250]; 

K_state=[8,9,10];
W=10;

if datatype==1
    if atlas==1
        N_node=100;
    else
        N_node=200;
    end
else
    N_node=100;
end
N_subj=100;


N_state=length(state_t);

group_adj=cell(N_subj,length(vari)); % group of adjacency matrix subject x DIIV

Q=cell(N_subj,length(vari));
K_g=cell(N_subj,length(vari)); % group of estimated K
z_g=cell(N_subj,length(vari)); % group of estimated z
true_latent_labels=cell(1,length(vari));

% K_chain=cell(N_subj,N_state);
% Z_chain=cell(N_subj,N_state);

nmi_score=cell(N_state,length(vari));

for ds=1:N_state
    for i=1:length(vari)
        nmi_score{ds,i}=zeros(N_subj,1);
    end
end
% K_chain=cell(N_subj,N_state);
% Z_chain=cell(N_subj,N_state);

for s=1:N_subj
    fprintf('Adjacency matrix of subject: %d\n',s)
    subid=num2str(subjects(s));
    for i=1:length(vari)
    [group_adj{s,i},true_latent,true_latent_labels{1,i}]=local_adj(datatype,atlas,subid,session_n,n_s,state_t,K_state,W,vari(i),hrf_ind);
    end
end

% -------------------------------------------------------------------------
% Adjacency matrix

M=10; % margin size
if atlas==1
    N=100; 
else
    N=200; 
end
group_adj_3D=cell(N_state,length(vari));
communities=zeros(N_node,N_subj,N_state);
Q_multi=zeros(1,N_state);


for ds=1:N_state
    for i=1:length(vari)
        group_adj_3D{ds,i}=cell(1,N_subj);
    end
end


for ds=1:N_state
    for i=1:length(vari)
        for s=1:N_subj
               group_adj_3D{ds,i}{1,s}=group_adj{s,i}{1,ds};
        end
    end
end

z_sub=cell(N_state,length(vari));
for ds=1:N_state
    for i=1:length(vari)
        z_sub{ds,i}=zeros(N_node,N_subj);
    end
end

for ds=1:N_state
    for i=1:length(vari)
        [z_sub{ds,i}] = multilayer_modularity(group_adj_3D{ds,i}(1,:), para{ds,i}(1,1), para{ds,i}(2,1));

    end
end


for ds=1:N_state
    for i=1:length(vari)
        
        for s=1:N_subj
            nmi_score{ds,i}(s,:)=nmi(z_sub{ds,i}(:,s),true_latent_labels{1,i}{s,1}(:,ds));
        end
        
    end
end


% Save the results of individual-level modelling

data_path = fileparts(mfilename('fullpath'));

group_path=fullfile(data_path,'Results/synthetic_Bayes_opt_MM/infer_opt');
save(group_path,'para','z_sub','true_latent_labels','nmi_score');











