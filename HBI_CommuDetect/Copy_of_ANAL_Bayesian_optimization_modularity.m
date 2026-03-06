% Bayesian optimization for modularity
% synthetic data analysis
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
Q=zeros(N_subj,N_state);
K_g=zeros(N_subj,N_state); % group of estimated K
z_g=cell(N_subj,N_state); % group of estimated z

for vari=10:10:30
    for s=1:N_subj
        fprintf('Adjacency matrix of subject: %d\n',s)
        subid=num2str(subjects(s));
        [group_adj{s,1},true_latent,true_latent_sub]=local_adj(datatype,atlas,subid,session_n,n_s,state_t,K_state,W,vari,hrf_ind);
    end
    
    % Hyperparameter optimization 
    for state=1:3
        fprintf('State: %d\n',state)
        opt_parameter(group_adj,true_latent_sub,n_s,N_subj,atlas,vari,state);
    end
   
end



function []=opt_parameter(group_adj,true_latent_sub,n_s,N_subj,atlas,vari,state)
  
    % Group-level Bayesian optimization (shared gamma)
    
    if atlas==1
        N=100; 
    else
        N=200; 
    end

    A_global = cell(N_subj,1);
    true_label = zeros(N,N_subj);
    K_true = zeros(N_subj,1);
    
    for s = 1:N_subj
        A_global{s} = group_adj{s,1}{1,state};
        true_label(:,s) = true_latent_sub{s,1}(:,state);
        K_true(s) = numel(unique(true_label(:,s)));
    end
    
    % pre-scan to tune lamda
    gamma_grid = linspace(1,2.5,50);
    
    NMIvals = zeros(size(gamma_grid));
    PenaltyVals = zeros(size(gamma_grid));
    
    for i = 1:length(gamma_grid)
    
        gamma = gamma_grid(i);
    
        NMI_sub = zeros(N_subj,1);
        penalty_sub = zeros(N_subj,1);
    
        for s = 1:N_subj
    
            z = modularity_und(A_global{s}, gamma);
    
            % NMI
            val = nmi(z, true_label(:,s));
            if isnan(val) || isinf(val)
                val = 0;
            end
            NMI_sub(s) = val;
    
            % penalty
            K_gamma = numel(unique(z));
            penalty_sub(s) = abs(K_gamma - K_true(s));
        end
    
        NMIvals(i) = mean(NMI_sub);
        PenaltyVals(i) = mean(penalty_sub);
    end
    
    range_NMI = max(NMIvals) - min(NMIvals);
    range_penalty = max(PenaltyVals) - min(PenaltyVals);
    
    if range_penalty < 1e-6
        lambda = 0;
    else
        lambda = 0.5 * range_NMI / range_penalty;
    end
    
    fprintf('Adaptive lambda (group-level) = %.4f\n', lambda);
    
    
    
    % -------------------------------------------------------------------------
    % Bayesian optimization
    
    objective = @(x) objective_fun_group(x, A_global, true_label, K_true, lambda);
    
    gammaVar = optimizableVariable('gamma',[1 2.5]);
    
    results = bayesopt(objective, gammaVar,...
        'MaxObjectiveEvaluations',30,...
        'IsObjectiveDeterministic',true,...
        'AcquisitionFunctionName','expected-improvement-plus');
    
    bestGamma = results.XAtMinObjective.gamma;
    
    fprintf('Best group-level gamma = %.4f\n', bestGamma);
    
    
    
    % -------------------------------------------------------------------------
    % extract data
    
    
    gammaTrace = results.XTrace.gamma;
    scoreTrace = -results.ObjectiveTrace;
    
    bestGamma  = results.XAtMinObjective.gamma;
    bestScore  = max(scoreTrace);
    
    bestSoFar  = cummax(scoreTrace);
    
    % ===============================
    % figures
    
    figure('Color','w','Position',[300 200 900 650])
    
    t = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    
    % Common style
    lw  = 1.8;
    ms  = 50;
    fs  = 11;
    
    % Gamma vs Objective
    nexttile
    scatter(gammaTrace, scoreTrace, ms, 'k','filled')
    hold on
    plot(bestGamma, bestScore,'p',...
        'MarkerSize',12,...
        'MarkerEdgeColor','r',...
        'MarkerFaceColor','r')
    
    xlabel('\gamma','FontSize',fs)
    ylabel('Objective Score','FontSize',fs)
    title('\gamma vs Objective','FontSize',fs+1)
    set(gca,'FontSize',fs,'LineWidth',1)
    box off
    grid on
    
    % Objective vs Iteration
    nexttile
    plot(scoreTrace,'-o',...
        'LineWidth',lw,...
        'MarkerSize',4)
    xlabel('Iteration','FontSize',fs)
    ylabel('Objective Score','FontSize',fs)
    title('Optimization Progress','FontSize',fs+1)
    set(gca,'FontSize',fs,'LineWidth',1)
    box off
    grid on
    
    % Convergence Curve
    nexttile
    plot(bestSoFar,'-','LineWidth',2.2)
    xlabel('Iteration','FontSize',fs)
    ylabel('Best Score So Far','FontSize',fs)
    title('Convergence','FontSize',fs+1)
    set(gca,'FontSize',fs,'LineWidth',1)
    box off
    grid on
    
    % Exploration Path
    nexttile
    plot(gammaTrace,'-o',...
        'LineWidth',lw,...
        'MarkerSize',4)
    xlabel('Iteration','FontSize',fs)
    ylabel('\gamma','FontSize',fs)
    title('Exploration Path','FontSize',fs+1)
    set(gca,'FontSize',fs,'LineWidth',1)
    box off
    grid on
    
    sgtitle('Bayesian Optimization Diagnostics','FontSize',14,'FontWeight','bold')
    
    data_path = fileparts(mfilename('fullpath'));
    
    NMI_indi_path=fullfile(data_path,['Results/synthetic_Bayes_opt_M','/','DIIV',num2str(vari),'/n',num2str(n_s),'/opt_para_',num2str(state)]);
    save(NMI_indi_path,'bestGamma')

    close all

end
% -------------------------------------------------------------------------
% nested function

function val = objective_fun_group(x, A_cell, true_label, K_true, lambda)

    gamma = x.gamma;

    N_subj = length(A_cell);
    score_sub = zeros(N_subj,1);

    for s = 1:N_subj

        z = modularity_und(A_cell{s}, gamma);

        % NMI
        nmi_score = nmi(z, true_label(:,s));
        if isnan(nmi_score) || isinf(nmi_score)
            nmi_score = 0;
        end

        % penalty
        K_gamma = numel(unique(z));
        penalty = abs(K_gamma - K_true(s));

        score_sub(s) = nmi_score - lambda * penalty;
    end

    val = -mean(score_sub);   % bayesopt minimizes
end







