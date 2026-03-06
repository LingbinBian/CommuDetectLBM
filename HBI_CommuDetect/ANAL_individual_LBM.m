% Community detection based on latent block model
% synthetic data analysis
%
% Version 1.0
% 27-Feb-2026
% Copyright (c) 2026, Lingbin Bian
% -------------------------------------------------------------------------
clear
clc
close all

ncore = str2double(getenv('SLURM_CPUS_PER_TASK'));
if isempty(gcp('nocreate'))
    parpool('local', ncore);
end
% -------------------------------------------------------------------------
% Load data
% Data type

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

hrf_ind=0;

state_t=[50,150,250]; 

K_state=[8,9,10];
W=10;

N_subj=100;
N=100; 
N_state=length(state_t);

atlas=1;

DIIV=[10,20,30];

opt_para_LBM=cell(3,3);
for state=1:3
    for n_diiv=1:length(DIIV)
        opt_para_LBM{state,n_diiv}=zeros(4,1);
    end
end

% load the optimal parameters
for state=1:3
    for n_diiv=1:length(DIIV)
    load(['/home/bianlb/Documents/CommuDetectLBM/HBI_CommuDetect/Results/synthetic_Bayes_opt_LBM','/','DIIV',num2str(DIIV(n_diiv)),'/n',num2str(n_s),'/','opt_para_',num2str(state)]);
    opt_para_LBM{state,n_diiv}(1,1)=bestK_max; %K_max
    opt_para_LBM{state,n_diiv}(2,1)=bestNu; %nu
    opt_para_LBM{state,n_diiv}(3,1)=bestRho; % rho
    opt_para_LBM{state,n_diiv}(4,1)=bestKappa_sq; % kappa_sq
    end
end

% K_max=23;
% nu=9.3596;
% rho=0.97034;
% kappa_sq=1.4192;

group_adj=cell(N_subj,length(DIIV)); % group of adjacency matrix, subject x DIIV

K_g=cell(N_subj,length(DIIV)); % group of estimated K
z_g=cell(N_subj,length(DIIV)); % group of estimated z

for s=1:N_subj
    for i=1:length(DIIV)
        z_g{s,i}=zeros(N,N_state);
    end
end

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
% Community detection

parfor s=1:N_subj
    fprintf('subject: %d\n',s)
    for i=1:3  
        for ds=1:N_state
          [z_g{s,i}(:,ds)] = CommuDetectLBM(group_adj{s,i}{1,ds},opt_para_LBM{ds,i}(1,1),opt_para_LBM{ds,i}(2,1),opt_para_LBM{ds,i}(3,1),opt_para_LBM{ds,i}(4,1));   % x,K_max,nu,rho,kappa_sq 
        end
    end
end

for ds=1:N_state
    for i=1:3       
        for s=1:N_subj
            nmi_score{ds,i}(s,:)=nmi(z_g{s,i}(:,ds),true_latent_labels{1,i}{s,1}(:,ds));
        end        
    end
end

% Save the results

group_path=fullfile('/home/bianlb/Documents/CommuDetectLBM/HBI_CommuDetect/Results/synthetic_Bayes_opt_LBM/infer_opt');
save(group_path,'opt_para_LBM','z_g','true_latent_labels','nmi_score');
   














