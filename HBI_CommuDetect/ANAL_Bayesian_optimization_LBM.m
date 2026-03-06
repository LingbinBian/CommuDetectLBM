% Bayesian optimization for latent block model (LBM) 
% synthetic data analysis
%
% Version 1.0
% 27-Feb-2026
% Copyright (c) 2026, Lingbin Bian
% -------------------------------------------------------------------------
clear
clc
close all

% parallel computing using SLURM high performance computer(HPC)
ncore = str2double(getenv('SLURM_CPUS_PER_TASK'));
if isempty(gcp('nocreate'))
    parpool('local', ncore);
end

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

atlas=1; % 1: ROI=100, 2: ROI=200
%vari=10; % variation of inter-individual variability
hrf_ind=0;
state_t=[50,150,250]; 
K_state=[8,9,10];
W=10;
N_subj=100;

group_adj=cell(N_subj,1); % group of adjacency matrix

for vari=10:10:30
    for s=1:N_subj
        fprintf('Adjacency matrix of subject: %d\n',s)
        subid=num2str(subjects(s));
        [group_adj{s,1},true_latent,true_latent_sub]=local_adj(datatype,atlas,subid,session_n,n_s,state_t,K_state,W,vari,hrf_ind);
    end
    
    
    for state=1:3
        fprintf('State: %d\n',state)
        opt_parameter(group_adj,true_latent_sub,n_s,N_subj,atlas,vari,state);
    end

end

% -------------------------------------------------------------------------
% nested function

function []=opt_parameter(group_adj,true_latent_sub,n_s,N_subj,atlas,vari,state)

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
    
    % -------------------------------------------------------------------------
    % Bayesian optimization
    % Joint optimization of nu, rho, kappa_sq, K_max
    % Objective: NMI 
    
      
    objective = @(x) objective_fun_group(x, A_global, true_label);
    
    vars = [
        optimizableVariable('nu',[2.1 15],'Transform','log')
        optimizableVariable('rho',[0.001 1],'Transform','log')
        optimizableVariable('kappa_sq',[0.1 20],'Transform','log')
        optimizableVariable('K_max',[6 25],'Type','integer')
        ];
    
    results = bayesopt(objective, vars,...
        'MaxObjectiveEvaluations',100,...
        'IsObjectiveDeterministic',false,...
        'AcquisitionFunctionName','expected-improvement-plus',...
        'Verbose',1);

    bestNu       = results.XAtMinObjective.nu;
    bestRho      = results.XAtMinObjective.rho;
    bestKappa_sq = results.XAtMinObjective.kappa_sq;
    bestK_max    = results.XAtMinObjective.K_max;
    bestScore    = results.MinObjective;

    fprintf('\nOptimal parameters (Group-level tuned):\n');
    fprintf('nu = %.4f\n', bestNu);
    fprintf('rho = %.4f\n', bestRho);
    fprintf('kappa_sq = %.4f\n', bestKappa_sq);
    fprintf('K_max = %d\n', bestK_max);
    fprintf('Best objective value = %.6f\n', bestScore);
    
    %data_path = fileparts(mfilename('fullpath'));
    
    para_path=fullfile(['/home/bianlb/Documents/CommuDetectLBM/HBI_CommuDetect/Results/synthetic_Bayes_opt_LBM','/','DIIV',num2str(vari),'/n',num2str(n_s),'/opt_para_',num2str(state)]);
    save(para_path,'bestNu','bestRho','bestKappa_sq','bestK_max','bestScore')
    
    
end
% ----------------------------------------------------------------------
% nested function

function val = objective_fun_group(x, A_cell, true_label)
% objective function

    nu       = x.nu;
    rho      = x.rho;
    kappa_sq = x.kappa_sq;
    K_max    = x.K_max;

    Nsub = length(A_cell);
    score_sub = zeros(Nsub,1);

    parfor s = 1:Nsub
        A_s = A_cell{s};  
        z_s = CommuDetectLBM(A_s,K_max,nu,rho,kappa_sq);
        score_sub(s) = nmi(z_s, true_label(:,s));
    end

    val = -mean(score_sub);
end





