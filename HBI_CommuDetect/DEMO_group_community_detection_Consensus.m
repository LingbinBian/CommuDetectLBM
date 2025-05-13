% Consensus clustering
% group-level analysis: modelling community memberships
%
% Version 1.0
% 20-April-2025
% Copyright (c) 2025, Lingbin Bian

clear
clc
close all

datatype=1;

if datatype==0
    
  % signal to noise ratio SNR=10log10(sigma_signal^2/sigma_noise^2);
  % sigma_signal=1, sigma_noise=n_s
  % n_s=0.3162;  % 10dB
  % n_s=0.5623;  % 5dB
  % n_s=1;       % 0dB
   n_s=1.7783;  % -5dB
  % n_s=3.1623;  % -10dB
  % n_s=5.6234;  % -15dB
  % n_s=10;     % -20dB
  % n_s=17.7828; % -25dB
    DIIV=10;
    load(['Results/synthetic_LBM/','DIIV',num2str(DIIV),'/','n',num2str(n_s),'/grouplevel_data.mat']);
elseif datatype==1
    session_n=2;
    if session_n==1
        load('Results/real/LR/grouplevel_data.mat');
    elseif session_n==2
        load('Results/real/RL/grouplevel_data.mat');
    end
end

% Modelling latent labels
mat_labels=zeros(N,N_subj); % temporary storage for a single state of latent labels

label_group_esti=zeros(N,N_state);
label_group_switched=zeros(N,N_state);
label_compare=cell(1,N_state);
vector_compare=zeros(N,2);

N_simu=100;
thre=0.5;
% Estimate assignment probability
for s=1:N_state
   for n=1:N_subj    
        mat_labels(:,n)=z_g{n,s};
   end
   K=max(max(mat_labels));
   label_group_esti(:,s)=CommuDetecConsensus(mat_labels,thre,N_simu);
   mat_labels=zeros(N,N_subj);
end


K_esti=max(label_group_esti);

if datatype==0
    for j=1:N_state
        vector_compare(:,2)=label_group_esti(:,j);
        vector_compare(:,1)=true_latent(:,j);
        vector_compare=labelswitch(vector_compare);
        label_compare{1,j}=vector_compare;        
    end
    for j=1:N_state
        label_group_switched(:,j)=label_compare{1,j}(:,2);
    end
end


% Visualize latent label vector

if datatype==0    
    for t=1:N_state
        %subplot(1,N_state,t)
        figure
        visual_labels(label_group_esti(:,t))
        title('Estimation','fontsize',16)
        set(gcf,'unit','normalized','position',[0.3,0.2,0.08,0.33]);
        saveas(gcf,['Results/synthetic_LBM/','DIIV',num2str(DIIV),'/n',num2str(n_s),'/Labels_esti_consensus_',num2str(t),'.fig'])

    end
     
    for t=1:N_state
        %subplot(1,N_state,t)
        figure
        visual_labels(label_compare{:,t}(:,2))
        title('Switched','fontsize',16)
        set(gcf,'unit','normalized','position',[0.3,0.2,0.08,0.33]);
        saveas(gcf,['Results/synthetic_LBM/','DIIV',num2str(DIIV),'/n',num2str(n_s),'/Labels_switched_consensus_',num2str(t),'.fig'])
    end

    for t=1:N_state
        %subplot(1,N_state,t)
        figure
        visual_labels(label_compare{:,t}(:,1))
        title('True','fontsize',16)
        set(gcf,'unit','normalized','position',[0.3,0.2,0.08,0.33]);
        % saveas(gcf,'Local_inference_synthetic/Labels_true.fig')
        saveas(gcf,['Results/synthetic_LBM/','DIIV',num2str(DIIV),'/n',num2str(n_s),'/Labels_true_consensus_',num2str(t),'.fig'])
    end
end

if datatype==1  
    for t=1:N_state
        %subplot(1,N_state,t)
        figure
        visual_labels(label_group_esti(:,t))
        title(['State',' ',num2str(t)],'fontsize',16) 
        set(gcf,'unit','normalized','position',[0.3,0.2,0.08,0.33]);
        if session_n==1
            saveas(gcf,['Results/real/LR/','Label_est_consensus_',num2str(t),'.fig'])
        elseif session_n==2
            saveas(gcf,['Results/real/RL/','Label_esti_consensus_',num2str(t),'.fig'])
        end

    end
end


data_path = fileparts(mfilename('fullpath'));
if datatype==0
   group_results_path=fullfile(data_path,['Results/synthetic_LBM/','DIIV',num2str(DIIV),'/n',num2str(n_s),'/grouplevel_results_',num2str(N_simu),'_consensus',]);
   save(group_results_path);
elseif datatype==1
   if session_n==1
       group_results_path=fullfile(data_path,['Results/real/LR/','grouplevel_results_',num2str(N_simu),'_consensus']);
   elseif session_n==2
       group_results_path=fullfile(data_path,['Results/real/RL/','grouplevel_results_',num2str(N_simu),'_consensus']);
   end
   save(group_results_path,'label_group_esti');    
end




