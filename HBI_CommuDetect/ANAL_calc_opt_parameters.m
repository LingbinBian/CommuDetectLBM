% Calculate the optimal parameters


% Version 1.0
% 16-Feb-2026
% Copyright (c) 2026, Lingbin Bian

clear
clc
close all

DIIV=[10,20,30];

n_s=0.3162;

N_simu=10;

opt_para_modularity=zeros(3,3); % state x DIIV
opt_para_multi_modularity=cell(3,3);
for state=1:3
    for n_diiv=1:length(DIIV)
        opt_para_multi_modularity{state,n_diiv}=zeros(2,1);

    end
end


% calculate the optimal parameters of modularity (synthetic data)
for state=1:3
    for n_diiv=1:length(DIIV)
        gamma=zeros(N_simu,1);
        for simu=1:N_simu
            load(['/Users/lingbinbian/Documents/MCD_revision_v2/CommuDetectLBM/HBI_CommuDetect/Results/synthetic_Bayes_opt_M/DIIV',num2str(DIIV(n_diiv)),'/n',num2str(n_s),'/opt_para_',num2str(state),'_simu_',num2str(simu),'.mat'])
            gamma(simu,1)=bestGamma;
        end
        opt_para_modularity(state,n_diiv)=mean(gamma);
    end
end

data_path = fileparts(mfilename('fullpath'));
para_path=fullfile(data_path,['Results/synthetic_Bayes_opt_M','/','opt_para_modularity']);
save(para_path,'opt_para_modularity')


% calculate the optimal parameters of multilayer modularity (synthetic data)
for state=1:3
    for n_diiv=1:length(DIIV)
        gamma=zeros(N_simu,1);
        omega=zeros(N_simu,1);
        for simu=1:N_simu
            load(['/Users/lingbinbian/Documents/MCD_revision_v2/CommuDetectLBM/HBI_CommuDetect/Results/synthetic_Bayes_opt_MM/DIIV',num2str(DIIV(n_diiv)),'/n',num2str(n_s),'/opt_para_',num2str(state),'_simu_',num2str(simu),'.mat'])
            gamma(simu,1)=bestGamma;
            omega(simu,1)=bestOmega;
        end
        opt_para_multi_modularity{state,n_diiv}(1,1)=mean(gamma);
        opt_para_multi_modularity{state,n_diiv}(2,1)=mean(omega);
    end
end

data_path = fileparts(mfilename('fullpath'));
para_path=fullfile(data_path,['Results/synthetic_Bayes_opt_MM','/','opt_para_multi_modularity']);
save(para_path,'opt_para_multi_modularity')



