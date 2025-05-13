% Community detection based on latent block model
% individual-level analysis
%
% Version 1.0
% 7-March-2024
% Copyright (c) 2024, Lingbin Bian
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
    atlas=2; % 1: Schaefer 100, 2: Kong 200
    if atlas==1
        N=100;
    else
        N=200;
    end
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
elseif datatype==0
    subjects=load('synthetic_id.txt');
    atlas='';
    session_n='';

  % signal to noise ratio SNR=10log10(sigma_signal^2/sigma_noise^2);
  % sigma_signal=1, sigma_noise=n_s
  % n_s=0.3162;  % 10dB
  % n_s=0.5623;  % 5dB
  % n_s=1;      
  % n_s=1.7783;  % -5dB
  % n_s=3.1623;  % -10dB
  % n_s=5.6234;  % -15dB
   n_s=10;      % -20dB
  % n_s=17.7828; % -25dB

    vari=20; % fixed true community structure
    hrf_ind=0;
    
    state_t=[50,150,250]; 

    K_state=[8,9,10];
    W=10;
end

N_subj=100;

N_state=length(state_t);

group_adj=cell(N_subj,1); % group of adjacency matrix
K_g=zeros(N_subj,N_state); % group of estimated K
z_g=cell(N_subj,N_state); % group of estimated z

K_chain=cell(N_subj,N_state);
z_chain=cell(N_subj,N_state);

for s=1:N_subj
    fprintf('Adjacency matrix of subject: %d\n',s)
    subid=num2str(subjects(s));
    [group_adj{s,1},true_latent,true_latent_sub]=local_adj(datatype,atlas,subid,session_n,n_s,state_t,K_state,W,vari,hrf_ind);
end

% -------------------------------------------------------------------------
% Adjacency matrix



for s=1:N_subj
    fprintf('Subject: %d\n----------------------------------------------------------\n',s)
    for ds=1:N_state
      fprintf('Subject: %d\n----------------------------------------------------------\n',s)
      fprintf('Condition: %d\n----------------------------------------------------------\n',ds)
      [z_g{s,ds},K_g(s,ds),z_chain{s,ds},K_chain{s,ds}] = CommuDetectLBM(group_adj{s,1}{1,ds});      
    end
end

M=10; % margin size

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
   group_path=fullfile(data_path,['Results/synthetic_LBM/','DIIV',num2str(vari),'/n',num2str(n_s),'/grouplevel_data']);
   save(group_path);
elseif datatype==1
    if atlas==1
       if session_n==1
           group_path=fullfile(data_path,['Results/real_LBM/LR/','grouplevel_data']);
       elseif session_n==2
           group_path=fullfile(data_path,['Results/real_LBM/RL/','grouplevel_data']);
       end
    else
       if session_n==1
           group_path=fullfile(data_path,['Results/real_Kong_LBM/LR/','grouplevel_data']);
       elseif session_n==2
           group_path=fullfile(data_path,['Results/real_Kong_LBM/RL/','grouplevel_data']);
       end
    end
   save(group_path);    
end

% Visualize latent label vectors ------------------------------------------
mat_labels=zeros(N,N_subj);
for s=1:N_state
    figure
    for n=1:N_subj    
        mat_labels(:,n)=z_g{n,s};  
    end
    imagesc(mat_labels(:,:))
    colormap(color_type(mat_labels(:,:)));
    colorbar_community_K(unique(mat_labels(:,:)))
    title(['State',' ',num2str(s)],'fontsize',16)
    xlabel('Subject','fontsize',16)
    ylabel('Node','fontsize',16)
    set(gca, 'linewidth', 1.2, 'fontsize', 16, 'fontname', 'times')
    mat_labels=zeros(N,N_subj);
    if datatype==0
        set(gcf,'unit','normalized','position',[0.3,0.2,0.15,0.35]);
    elseif datatype==1
        set(gcf,'unit','normalized','position',[0.3,0.2,0.15,0.35]);
    end
    if datatype==0
        saveas(gcf,['Results/synthetic_LBM/','DIIV',num2str(vari),'/n',num2str(n_s),'/group_latent_labels_',num2str(s),'.fig'])
    elseif datatype==1
        if atlas==1
            if session_n==1
                saveas(gcf,['Results/real_LBM/LR/','group_latent_labels_',num2str(s),'.fig'])
            elseif session_n==2
                saveas(gcf,['Results/real_LBM/RL/','group_latent_labels_',num2str(s),'.fig'])
            end
        else
            if session_n==1
                saveas(gcf,['Results/real_Kong_LBM/LR/','group_latent_labels_',num2str(s),'.fig'])
            elseif session_n==2
                saveas(gcf,['Results/real_Kong_LBM/RL/','group_latent_labels_',num2str(s),'.fig'])
            end
        end
    end
end










