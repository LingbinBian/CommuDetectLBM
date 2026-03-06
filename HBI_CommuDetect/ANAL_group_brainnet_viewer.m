% group-level brain network visualization
% Write memberships (labels) of nodes, and edges (adjacency matrix) to .node and .edge files for BrainNet Viewer.

% Version 1.0
% 7-Jun-2020
% Copyright (c) 2020, Lingbin Bian
clear
clc

atlas=2;
L_state=11;

session_n=1;

if atlas==1
    if session_n==1
       cd 'Results/real_LBM/LR'
    end
    if session_n==2
       cd 'Results/real_LBM/RL'
    end
else
    if session_n==1
       cd 'Results/real_Kong_LBM/LR'
    end
    if session_n==2
       cd 'Results/real_Kong_LBM/RL'
    end
end
load('grouplevel_data.mat');  
load('grouplevel_results.mat');
spar_level=5;
if atlas==1
    N=100;
else
    N=200;
end
% esti_grouplabel=zeros(N,11);
% for t=1:L_state
%     for i=1:N
%         for j=1:12
%             if R_esti_max{1,t}(i,j)~=0
%                 esti_grouplabel(i,t)=j;
%             end
%         end
%     end      
% end


% data_path = fileparts(mfilename('fullpath'));
% if datatype==0
%    group_path=fullfile(data_path,['Results/synthetic/','n',num2str(n_s),'/grouplevel_est']);
%    save(group_path);
% elseif datatype==1
%    if session_n==1
%        group_path=fullfile(data_path,['Results/real/LR/','grouplevel_est']);
%    elseif session_n==2
%        group_path=fullfile(data_path,['Results/real/RL/','grouplevel_est']);
%    end
%    save(group_path);    
% end

if datatype==1
    if session_n==1
        load('LR_meanconnectivity.mat')
    elseif session_n==2
        load('RL_meanconnectivity.mat')
    end
end

for t=1:3
    Node_ROI=dlmread(['sparsity_',num2str(spar_level),'/Node_ROI_t',num2str(state_t(t)),'.node']);
    Node_ROI(:,4)=label_group_esti(:,t);
    if atlas==2
        Node_ROI(:,6)=1:200;
    end
    dlmwrite(['sparsity_',num2str(spar_level),'/Node_ROI_t',num2str(state_t(t)),'.node'],Node_ROI,'delimiter','\t')
    dlmwrite(['sparsity_',num2str(spar_level),'/Weighted_t',num2str(state_t(t)),'.edge'],ave_adj{t},'delimiter','\t')
   
     BrainNet_MapCfg(['sparsity_',num2str(spar_level),'/Node_ROI_t',num2str(state_t(t)),'.node'],...
         ['sparsity_',num2str(spar_level),'/Weighted_t',num2str(state_t(t)),'.edge'],...
         '/Users/lingbinbian/Documents/PhD_Project/BrainNet_viewer/BrainNetViewer_20191031/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv',...
         ['sparsity_',num2str(spar_level),'/Brainnet_setup_t',num2str(state_t(t)),'.mat'],...
         ['sparsity_',num2str(spar_level),'/network_t',num2str(state_t(t)),'.fig']);
     saveas(gcf,['sparsity_',num2str(spar_level),'/network_t',num2str(state_t(t)),'.fig'])
     saveas(gcf,['sparsity_',num2str(spar_level),'/network_t',num2str(state_t(t)),'.svg'])
   
end
cd ..

% module colors in brainnet viewer
%[0.69,0.19,0.38; 1.00,0.89,0.52; 0.25,0.41,0.88; 0.00,0.79,0.34; 0.63,0.32,0.18; 0.50,0.54,0.53; 0.53,0.81,0.92; 0.63,0.40,0.83; 0.74,0.56,0.56;0.16,0.14,0.13;1.00,0.50,0.31;1.00,0.00,0.00;0.39,0.19,0.38;0.50,0.89,0.52;0.12,0.41,0.88;0.31,0.32,0.18]