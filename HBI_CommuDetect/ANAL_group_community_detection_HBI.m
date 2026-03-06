% Hierarchical Bayesian inference
% group-level analysis: modelling community memberships
%
% Version 1.0
% 4-May-2025
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
  % n_s=1.7783;  % -5dB
  % n_s=3.1623;  % -10dB
  % n_s=5.6234;  % -15dB
  %  n_s=10;     % -20dB
  % n_s=17.7828; % -25dB
    DIIV=10;
    load(['Results/synthetic_LBM/','DIIV',num2str(DIIV),'/','n',num2str(n_s),'/grouplevel_data.mat']);
elseif datatype==1
    session_n=1;
    atlas=2;
    if atlas==1
        if session_n==1
            load('Results/real_LBM/LR/grouplevel_data.mat');
        elseif session_n==2
            load('Results/real_LBM/RL/grouplevel_data.mat');
        end
    else
        if session_n==1
            load('Results/real_Kong_LBM/LR/grouplevel_data.mat');
        elseif session_n==2
            load('Results/real_Kong_LBM/RL/grouplevel_data.mat');
        end

    end
end

% Modelling latent labels
mat_labels=zeros(N,N_subj); % temporary storage for a single state of latent labels
R_esti=cell(1,N_state); % label assignment probability matrix
R_esti_max=cell(1,N_state); % maximum label assignment probability matrix
label_group_esti=zeros(N,N_state);
label_group_switched=zeros(N,N_state);
label_compare=cell(1,N_state);
vector_compare=zeros(N,2);

N_simu=1000;

% Estimate assignment probability
for s=1:N_state
   for n=1:N_subj    
        mat_labels(:,n)=z_g{n,s};
   end
   K=max(max(mat_labels));
   R_esti{1,s}=assign_esti(mat_labels,K,N_simu);
   mat_labels=zeros(N,N_subj);
end

for s=1:N_state
   R_esti_max{s}=mat_maxrow(R_esti{s});
   R_esti_max{s}(:,all(R_esti_max{s}==0,1))=[];
end

for s=1:N_state
    [N,K_assign]=size(R_esti_max{s});
    for i=1:N
      for j=1:K_assign
          if R_esti_max{s}(i,j)~=0  
             label_group_esti(i,s)=j;
          end
      end
    end
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


% Visualize assignment probability

for s=1:N_state
    figure
    imagesc(R_esti{s})
    colormap(sky);
    colorbar
    title(['Condition',' ',num2str(s)],'fontsize',16) 
    xlabel('k','fontsize',16) 
    ylabel('Node','fontsize',16)      
    set(gca, 'linewidth', 1.2, 'fontsize', 16, 'fontname', 'times')
    if datatype==0
        set(gcf,'unit','normalized','position',[0.3,0.2,0.15,0.35]);
    elseif datatype==1
        set(gcf,'unit','normalized','position',[0.3,0.2,0.15,0.35]);
    end
    
    if datatype==0
        saveas(gcf,['Results/synthetic_LBM/','DIIV',num2str(DIIV),'/n',num2str(n_s),'/assign_prob_state_',num2str(s),'.fig'])
    elseif datatype==1
        if atlas==1
            if session_n==1
                saveas(gcf,['Results/real_LBM/LR/','assign_prob_state_',num2str(s),'.fig'])
            elseif session_n==2
                saveas(gcf,['Results/real_LBM/RL/','assign_prob_state_',num2str(s),'.fig'])
            end
        else
            if session_n==1
                saveas(gcf,['Results/real_Kong_LBM/LR/','assign_prob_state_',num2str(s),'.fig'])
            elseif session_n==2
                saveas(gcf,['Results/real_Kong_LBM/RL/','assign_prob_state_',num2str(s),'.fig'])
            end
        end
    end
end

% Visualize assignment probability (Max)

for s=1:N_state
   if datatype==0  
      figure      
   elseif datatype==1
      figure
   end 
   
   imagesc(R_esti_max{s})
   colormap(sky);
   colorbar
   title(['Condition',' ',num2str(s)],'fontsize',16) 
   xlabel('k','fontsize',16) 
   ylabel('Node','fontsize',16)    
   set(gca, 'linewidth', 1.2, 'fontsize', 16, 'fontname', 'times')


   if datatype==0
    set(gcf,'unit','normalized','position',[0.3,0.2,0.15,0.35]);
   elseif datatype==1
    set(gcf,'unit','normalized','position',[0.3,0.2,0.15,0.35]);
   end
%    
   if datatype==0
       saveas(gcf,['Results/synthetic_LBM/','DIIV',num2str(DIIV),'/n',num2str(n_s),'/assign_prob_max_state_',num2str(s),'.fig'])
   elseif datatype==1
       if atlas==1
           if session_n==1
               saveas(gcf,['Results/real_LBM/LR/','assign_prob_max_state_',num2str(s),'.fig'])
           elseif session_n==2
               saveas(gcf,['Results/real_LBM/RL/','assign_prob_max_state_',num2str(s),'.fig'])
           end
       else
           if session_n==1
               saveas(gcf,['Results/real_Kong_LBM/LR/','assign_prob_max_state_',num2str(s),'.fig'])
           elseif session_n==2
               saveas(gcf,['Results/real_Kong_LBM/RL/','assign_prob_max_state_',num2str(s),'.fig'])
           end

       end
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
        saveas(gcf,['Results/synthetic_LBM/','DIIV',num2str(DIIV),'/n',num2str(n_s),'/Labels_esti_',num2str(t),'.fig'])

    end
     
    for t=1:N_state
        %subplot(1,N_state,t)
        figure
        visual_labels(label_compare{:,t}(:,2))
        title('Switched','fontsize',16)
        set(gcf,'unit','normalized','position',[0.3,0.2,0.08,0.33]);
        saveas(gcf,['Results/synthetic_LBM/','DIIV',num2str(DIIV),'/n',num2str(n_s),'/Labels_switched_',num2str(t),'.fig'])
    end

    for t=1:N_state
        %subplot(1,N_state,t)
        figure
        visual_labels(label_compare{:,t}(:,1))
        title('True','fontsize',16)
        set(gcf,'unit','normalized','position',[0.3,0.2,0.08,0.33]);
        % saveas(gcf,'Local_inference_synthetic/Labels_true.fig')
        saveas(gcf,['Results/synthetic_LBM/','DIIV',num2str(DIIV),'/n',num2str(n_s),'/Labels_true_',num2str(t),'.fig'])
    end
end

if datatype==1  
    for t=1:N_state
        %subplot(1,N_state,t)
        figure
        visual_labels(label_group_esti(:,t))
        title(['State',' ',num2str(t)],'fontsize',16) 
        set(gcf,'unit','normalized','position',[0.3,0.2,0.08,0.33]);
        if atlas==1
            if session_n==1
                saveas(gcf,['Results/real_LBM/LR/','Label_est_',num2str(t),'.fig'])
            elseif session_n==2
                saveas(gcf,['Results/real_LBM/RL/','Label_esti_',num2str(t),'.fig'])
            end
        else
            if session_n==1
                saveas(gcf,['Results/real_Kong_LBM/LR/','Label_est_',num2str(t),'.fig'])
            elseif session_n==2
                saveas(gcf,['Results/real_Kong_LBM/RL/','Label_esti_',num2str(t),'.fig'])
            end
        end

    end
end


data_path = fileparts(mfilename('fullpath'));
if datatype==0
   group_results_path=fullfile(data_path,['Results/synthetic_LBM/','DIIV',num2str(DIIV),'/n',num2str(n_s),'/grouplevel_results']);
   save(group_results_path);
elseif datatype==1
   if atlas==1 
       if session_n==1
           group_results_path=fullfile(data_path,['Results/real_LBM/LR/','grouplevel_results']);
       elseif session_n==2
           group_results_path=fullfile(data_path,['Results/real_LBM/RL/','grouplevel_results']);
       end
   else
       if session_n==1
           group_results_path=fullfile(data_path,['Results/real_Kong_LBM/LR/','grouplevel_results']);
       elseif session_n==2
           group_results_path=fullfile(data_path,['Results/real_Kong_LBM/RL/','grouplevel_results']);
       end
   end
   save(group_results_path,'R_esti_max','label_group_esti');    
end




