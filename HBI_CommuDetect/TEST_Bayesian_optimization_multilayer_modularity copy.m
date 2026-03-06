% Bayesian optimization for multilayer modularity
%
%
% Version 1.0
% 16-Feb-2026
% Copyright (c) 2026, Lingbin Bian
% -------------------------------------------------------------------------
clear
clc
close all
% -------------------------------------------------------------------------
% Load data
% Data type
datatype=0;   % 1: real data, 0: synthetic data

if datatype==1
    subjects=load('subject.txt');
    atlas=1;
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
    % gamma=2;
    atlas=1;
    vari=10; % fixed true community structure
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

for state=3:3
    fprintf('State: %d\n',state)
    opt_parameter(group_adj,true_latent_sub,n_s,N_subj,N_state,atlas,vari,state);
end


function []=opt_parameter(group_adj,true_latent_sub,n_s,N_subj,N_state,atlas,vari,state)

    M=10; % margin size
    if atlas==1
        N_node=100; 
    else
        N_node=200; 
    end
    
    group_adj_3D=cell(N_state,N_subj);
    communities=zeros(N_node,N_subj,N_state);
    Q_multi=zeros(1,N_state);
    
    for s=1:N_subj
        fprintf('subject: %d\n',s)
        for ds=1:N_state
           group_adj_3D{ds,s}=group_adj{s,1}{1,ds};
        end
    end
    
    % -------------------------------------------------------------------------
    % Bayesian optimization
    % Joint optimization of gamma and omega
    % Objective: NMI - lambda * K penalty
    
    
    A_global=group_adj_3D(state,:);    % multilayer adjacency
    true_label=zeros(N_node,N_subj);
    K=zeros(1,N_subj);
    
    for s=1:N_subj
        true_label(:,s)=true_latent_sub{s,1}(:,state);
        K(s)=numel(unique(true_label(:,s)));
    end
    
    % ----------------------------------------------------------
    % ----------------------------------------------------------
    % Pre-scan (gamma, omega) to estimate adaptive lambda
    % Stable & reviewer-safe version
    % ----------------------------------------------------------
    
    % gamma_grid = linspace(1,2.5,30);
    % omega_grid = linspace(0,1,30);
    % 
    % NMIvals = zeros(length(gamma_grid), length(omega_grid));
    % PenaltyVals = zeros(length(gamma_grid), length(omega_grid));
    % 
    % N_subj = size(true_label,2);
    % 
    % for i = 1:length(gamma_grid)
    %     for j = 1:length(omega_grid)
    % 
    %         z = multilayer_modularity(A_global,...
    %                                   gamma_grid(i),...
    %                                   omega_grid(j));
    % 
    %         NMI_sub = zeros(N_subj,1);
    %         penalty_sub = zeros(N_subj,1);
    % 
    %         for s = 1:N_subj
    % 
    %             z_s = z(:,s);
    %             true_s = true_label(:,s);
    % 
    %             % ----- NMI -----
    %             NMI_sub(s) = nmi(z_s, true_s);
    % 
    %             % ----- normalized K penalty -----
    %             K_para = numel(unique(z_s));
    %             K_true = numel(unique(true_s));
    % 
    %             penalty_sub(s) = abs(K_para - K_true) / (K_true + eps);
    %         end
    % 
    %         NMIvals(i,j) = mean(NMI_sub);
    %         PenaltyVals(i,j) = mean(penalty_sub);
    %     end
    % end
    % 
    % % ----------------------------------------------------------
    % % Robust lambda estimation (std-based scaling)
    % % ----------------------------------------------------------
    % 
    % std_NMI = std(NMIvals(:));
    % std_penalty = std(PenaltyVals(:));
    % 
    % if std_penalty < 1e-6
    %     lambda = 0;
    % else
    %     lambda = std_NMI / std_penalty;
    % end
    % 
    % fprintf('Adaptive lambda (robust) = %.6f\n', lambda);

    
    % ----------------------------------------------------------
    % Bayesian optimization    
    objective = @(x) objective_fun_group(x, A_global, true_label);    
    vars = [
        optimizableVariable('gamma',[1 2.5])
        optimizableVariable('omega',[0 1])
        ];   
    results = bayesopt(objective, vars,...
        'MaxObjectiveEvaluations',30,...
        'IsObjectiveDeterministic',true,...
        'AcquisitionFunctionName','expected-improvement-plus');   
    bestGamma = results.XAtMinObjective.gamma;
    bestOmega = results.XAtMinObjective.omega;
    
    fprintf('\nOptimal parameters (Group-level tuned):\n');
    fprintf('gamma = %.4f\n', bestGamma);
    fprintf('omega = %.4f\n', bestOmega);
    
    data_path = fileparts(mfilename('fullpath'));
    
    NMI_indi_path=fullfile(data_path,['Results/synthetic_Bayes_opt_MM','/','DIIV',num2str(vari),'/n',num2str(n_s),'/opt_para_',num2str(state)]);
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
end
% ----------------------------------------------------------------------
% Objective Function
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

        % K penalty
        %K_para = numel(unique(z_s));
        %K_true = numel(unique(true_s));

        %penalty = abs(K_para - K_true);

        % normalized penalty
        %penalty_norm = penalty / (K_true + eps);

        %score_sub(s) = nmi_score - lambda * penalty_norm;
        score_sub(s) = nmi_score;
    end

    val = -mean(score_sub);
end



