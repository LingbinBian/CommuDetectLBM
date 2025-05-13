% Community detection using modularity
% individual-level analysis
%
% Version 1.0
% 6-March-2025
% Copyright (c) 2025, Lingbin Bian
% -------------------------------------------------------------------------
clear
clc
close all
% -------------------------------------------------------------------------
% Load data
% Data type
datatype=1;   % 1: real data, 0: synthetic data

if datatype==1
    subjects=load('subject.txt');
    atlas=2;
    session_n=1;
    n_s='';
    vari='';
    hrf_ind='';
    if session_n==1
        state_t=[41,76,107,140,175,206,239,278,306,334,375];
    elseif session_n==2
        state_t=[49,77,99,139,175,209,236,275,305,334,376];
    end
    K_state=[6,6,6,6,6,6,6,6,6,6,6];  % make an assumption of K
    W=10;
    gamma=2;
elseif datatype==0
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

   % modularity resolution gamma=1:0.1:2
    gamma=2;
    
    vari=30; % fixed true community structure
    hrf_ind=0;
    
    state_t=[50,150,250]; 

    K_state=[8,9,10];
    W=10;
end

N_subj=100;

N_state=length(state_t);

group_adj=cell(N_subj,1); % group of adjacency matrix
Q=zeros(N_subj,N_state);
K_g=zeros(N_subj,N_state); % group of estimated K
z_g=cell(N_subj,N_state); % group of estimated z

% K_chain=cell(N_subj,N_state);
% Z_chain=cell(N_subj,N_state);

for s=1:N_subj
    fprintf('Adjacency matrix of subject: %d\n',s)
    subid=num2str(subjects(s));
    [group_adj{s,1},true_latent,true_latent_sub]=local_adj(datatype,atlas,subid,session_n,n_s,state_t,K_state,W,vari,hrf_ind);
end

% -------------------------------------------------------------------------
% Adjacency matrix

M=10; % margin size
if atlas==1
    N=100; 
else
    N=200; 
end

for s=1:N_subj
    fprintf('subject: %d\n',s)
    for ds=1:N_state
      [z_g{s,ds},Q(s,ds)] = modularity_und(group_adj{s,1}{1,ds},gamma);   
      K_g(s,ds)=max(z_g{s,ds});
    end
end

% Label switching
group_labels=zeros(N,N_subj*N_state);
for ds=1:N_state
    for s=1:N_subj
       group_labels(:,N_subj*(ds-1)+s)=z_g{s,ds};
    end
end
group_labels=labelswitch(group_labels);
for ds=1:N_state
    for s=1:N_subj
       z_g{s,ds}=group_labels(:,N_subj*(ds-1)+s);   
    end
end

% Save the results of individual-level modelling

data_path = fileparts(mfilename('fullpath'));
if datatype==0

   group_path=fullfile(data_path,['Results/synthetic_modularity/','DIIV',num2str(vari),'/n',num2str(n_s),'/',num2str(gamma),'/grouplevel_data']);
   save(group_path);
   
elseif datatype==1
   if atlas==1
       if session_n==1
           group_path=fullfile(data_path,['Results/real_modularity/',num2str(gamma),'/LR/grouplevel_data']);
       elseif session_n==2
           group_path=fullfile(data_path,['Results/real_modularity/',num2str(gamma),'/RL/grouplevel_data']);
       end
       save(group_path);   
   else
       if session_n==1
           group_path=fullfile(data_path,['Results/real_Kong_modularity/',num2str(gamma),'/LR/grouplevel_data']);
       elseif session_n==2
           group_path=fullfile(data_path,['Results/real_Kong_modularity/',num2str(gamma),'/RL/grouplevel_data']);
       end
       save(group_path); 
   end
end












