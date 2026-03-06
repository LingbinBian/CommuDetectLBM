% Bayesian optimization for multilayer modularity
% synthetic data analysis
%
% Version 1.0
% 16-Feb-2026
% Copyright (c) 2026, Lingbin Bian
% -------------------------------------------------------------------------
addpath('/home/bianlb/Documents/CommuDetectLBM/HBI_CommuDetect/GenLouvain-master')

clear
clc

% -------------------------------------------------------------------------
% Load data
% Data type
datatype=0;   % 1: real data, 0: synthetic data


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

hrf_ind=0;

state_t=[50,150,250]; 
K_state=[8,9,10];
W=10;

N_subj=100;

N_state=length(state_t);

group_adj=cell(N_subj,1); % group of adjacency matrix
z_g=cell(N_subj,N_state); % group of estimated z

for vari=10:10:30
    for s=1:N_subj
        fprintf('Adjacency matrix of subject: %d\n',s)
        subid=num2str(subjects(s));
        [group_adj{s,1},true_latent,true_latent_sub]=local_adj(datatype,atlas,subid,session_n,n_s,state_t,K_state,W,vari,hrf_ind);
    end
    
    % Hyperparameter optimization 
    for state=3:3
        fprintf('State: %d\n',state)
        opt_parameter(group_adj,true_latent_sub,n_s,N_subj,N_state,atlas,vari,state);
    end

end

function []=opt_parameter(group_adj,true_latent_sub,n_s,N_subj,N_state,atlas,vari,state)

    if atlas==1
        N_node=100; 
    else
        N_node=200; 
    end
    
    group_adj_3D=cell(N_state,N_subj);
    
    for s=1:N_subj
        fprintf('subject: %d\n',s)
        for ds=1:N_state
           group_adj_3D{ds,s}=group_adj{s,1}{1,ds};
        end
    end

    % Bayesian optimization
    % Joint optimization of gamma and omega
    % Objective: NMI 
    
    A_global=group_adj_3D(state,:);    % multilayer adjacency
    true_label=zeros(N_node,N_subj);

    
    for s=1:N_subj
        true_label(:,s)=true_latent_sub{s,1}(:,state);
       
    end
    
    % ----------------------------------------------------------
      
    objective = @(x) objective_fun_group(x, A_global, true_label);    
    vars = [
        optimizableVariable('gamma',[1 2.5])
        optimizableVariable('omega',[0 0.1])
        ];   
    results = bayesopt(objective, vars,...
        'MaxObjectiveEvaluations',100,...
        'IsObjectiveDeterministic',true,...
        'AcquisitionFunctionName','expected-improvement-plus');   
    bestGamma = results.XAtMinObjective.gamma;
    bestOmega = results.XAtMinObjective.omega;
    
    fprintf('\nOptimal parameters (Group-level tuned):\n');
    fprintf('gamma = %.4f\n', bestGamma);
    fprintf('omega = %.4f\n', bestOmega);
    
    
    NMI_indi_path=fullfile(['/home/bianlb/Documents/CommuDetectLBM/HBI_CommuDetect/Results/synthetic_Bayes_opt_MM','/','DIIV',num2str(vari),'/n',num2str(n_s),'/opt_para_',num2str(state)]);
    save(NMI_indi_path,'bestGamma','bestOmega')
    
    
    % ----------------------------------------------------------
    % Convergence Plot    
    scoreTrace = -results.ObjectiveTrace;
    bestSoFar  = cummax(scoreTrace);    
    figure('Color','w')
    plot(bestSoFar,'LineWidth',2)
    xlabel('Iteration')
    ylabel('Best Objective So Far')
    title('Joint Optimization Convergence')
    grid on
          
    % ----------------------------------------------------------
    % Visualization    
    scoreTrace = -results.ObjectiveTrace;
    bestSoFar  = cummax(scoreTrace);  
    figure('Color','w')
    plot(bestSoFar,'LineWidth',2)
    xlabel('Iteration')
    ylabel('Best Mean Score So Far')
    title('Group-Level Optimization Convergence')
    grid on
    close all
end
% ----------------------------------------------------------------------
% Objective function
function val = objective_fun_group(x, A, true_label)

    gamma = x.gamma;
    omega = x.omega;

    z = multilayer_modularity(A, gamma, omega);

    Nsub = size(z,2);
    score_sub = zeros(Nsub,1);

    for s = 1:Nsub
        z_s   = z(:,s);
        true_s = true_label(:,s);
        % NMI
        nmi_score = nmi(z_s, true_s);          
        score_sub(s) = nmi_score;
    end
    val = -mean(score_sub);
end



