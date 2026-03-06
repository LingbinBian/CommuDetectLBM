% Calculate the optimal parameters (real fMRI data)

% Version 1.0
% 23-Feb-2026
% Copyright (c) 2026, Lingbin Bian

clear
clc
close all

N_simu=10;
% -------------------------------------------------------------------------
% modularity
% Schaefer's atlas
opt_para_M_indi=zeros(N_simu,6);

for c=1:6  % combination of working memory paradigm blocks
    for simu=1:N_simu
        load(['/Users/lingbinbian/Documents/MCD_revision_v2/CommuDetectLBM/HBI_CommuDetect/Results/real_Bayes_opt_M/combination',num2str(c),'_simu_',num2str(simu),'.mat'])
        opt_para_M_indi(simu,c)=bestGamma;    
    end
end

opt_para_M=mean(opt_para_M_indi);

data_path = fileparts(mfilename('fullpath'));
para_path=fullfile(data_path,['Results/real_Bayes_opt_M','/','opt_para_modularity']);
save(para_path,'opt_para_M')

% Kong's atlas
opt_para_M_indi_Kong=zeros(N_simu,6);

for c=1:6  % combination of working memory paradigm blocks
    for simu=1:N_simu
        load(['/Users/lingbinbian/Documents/MCD_revision_v2/CommuDetectLBM/HBI_CommuDetect/Results/real_Kong_Bayes_opt_M/combination',num2str(c),'_simu_',num2str(simu),'.mat'])
        opt_para_M_indi_Kong(simu,c)=bestGamma;    
    end
end

opt_para_M_Kong=mean(opt_para_M_indi_Kong);

data_path = fileparts(mfilename('fullpath'));
para_path=fullfile(data_path,['Results/real_Kong_Bayes_opt_M','/','opt_para_modularity']);
save(para_path,'opt_para_M_Kong')
% -------------------------------------------------------------------------
% multilayer modularity
% Schaefer's atlas
opt_para_MM=cell(2,6);

for c=1:6  % combination
    gamma=zeros(N_simu,1);
    omega=zeros(N_simu,1);
    for simu=1:N_simu

        load(['/Users/lingbinbian/Documents/MCD_revision_v2/CommuDetectLBM/HBI_CommuDetect/Results/real_Bayes_opt_MM/combination',num2str(c),'_simu_',num2str(simu),'.mat'])
        gamma(simu,1)=bestGamma; 
        omega(simu,1)=bestOmega; 
    end
    opt_para_MM{1,c}=mean(gamma);
    opt_para_MM{2,c}=mean(omega);
end

data_path = fileparts(mfilename('fullpath'));
para_path=fullfile(data_path,['Results/real_Bayes_opt_MM','/','opt_para_multi_modularity']);
save(para_path,'opt_para_MM')

% Kong's atlas

opt_para_MM_Kong=cell(2,6);

for c=1:6  % combination
    gamma=zeros(N_simu,1);
    omega=zeros(N_simu,1);
    for simu=1:N_simu

        load(['/Users/lingbinbian/Documents/MCD_revision_v2/CommuDetectLBM/HBI_CommuDetect/Results/real_Kong_Bayes_opt_MM/combination',num2str(c),'_simu_',num2str(simu),'.mat'])
        gamma(simu,1)=bestGamma; 
        omega(simu,1)=bestOmega; 
    end
    opt_para_MM_Kong{1,c}=mean(gamma);
    opt_para_MM_Kong{2,c}=mean(omega);
end

data_path = fileparts(mfilename('fullpath'));
para_path=fullfile(data_path,['Results/real_Kong_Bayes_opt_MM','/','opt_para_multi_modularity']);
save(para_path,'opt_para_MM_Kong')






