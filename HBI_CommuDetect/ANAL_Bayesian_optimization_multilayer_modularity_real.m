% Bayesian optimization for multilayer modularity (real data)
%
% Version 1.0
% 21-Feb-2026
% Copyright (c) 2026, Lingbin Bian
% -------------------------------------------------------------------------
clear
clc
close all
% -------------------------------------------------------------------------
% Load data
% Data type

datatype=1;   % 1: real data, 0: synthetic data

subjects=load('subject.txt');
atlas=1; % 1: Schaefer 100, 2: Kong 200
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

N_subj=100;

N_state=length(state_t);

group_adj=cell(N_subj,1); % group of adjacency matrix
Q=zeros(N_subj,N_state);
K_g=zeros(N_subj,N_state); % group of estimated K
z_g=cell(N_subj,N_state); % group of estimated z


for s=1:N_subj
    fprintf('Adjacency matrix of subject: %d\n',s)
    subid=num2str(subjects(s));
    [group_adj{s,1},true_latent,true_latent_sub]=local_adj(datatype,atlas,subid,session_n,n_s,state_t,K_state,W,vari,hrf_ind);
end


%--------------------------------------------------------------------------
% 1:2-back tool, 2:0-back body, 3:fixation, 4:2-back face, 5:0-back tool, 6:fixation, 7:2-back body,
% 8:2-back place, 9:fixation, 10:0-back face, 11:0-back place, 12:fixation

% 2-back    
for c_comb=1:6
    switch c_comb
        case 1
            c=[1,4,7,8];       
        case 2  
            c=[1,7,4,8];           
        case 3
            c=[1,8,4,7];   
        case 4
            c=[2,5,10,11];
        case 5
            c=[2,10,5,11];
        case 6
            c=[2,11,5,10];
            
    end
    A_global = cell(N_subj,4);
    for s = 1:N_subj
        A_global{s,1} = group_adj{s,1}{1,c(1)};
        A_global{s,2} = group_adj{s,1}{1,c(2)};
        A_global{s,3} = group_adj{s,1}{1,c(3)};
        A_global{s,4} = group_adj{s,1}{1,c(4)};     
    end                               
    opt_parameter(A_global,atlas,c_comb);
  
end

function []=opt_parameter(A_global,atlas,c_comb)
    % Group-level Bayesian optimization 
        
    objective = @(x) objective_fun_group(x, A_global);    
    vars = [
        optimizableVariable('gamma',[2 4])
        optimizableVariable('omega',[0 0.05])
        ];   
    results = bayesopt(objective, vars,...
        'MaxObjectiveEvaluations',100,...
        'IsObjectiveDeterministic',false,...
        'AcquisitionFunctionName','expected-improvement-plus');   
    bestGamma = results.XAtMinObjective.gamma;
    bestOmega = results.XAtMinObjective.omega;
    bestScore = results.MinObjective;
   
    fprintf('\nOptimal parameters (Group-level tuned):\n');
    fprintf('gamma = %.4f\n', bestGamma);
    fprintf('omega = %.4f\n', bestOmega);
    fprintf('Best objective value = %.6f\n', bestScore);
    
    
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

    if atlas==1
       
        NMI_indi_path=fullfile(['/home/bianlb/Documents/CommuDetectLBM/HBI_CommuDetect/Results/real_Bayes_opt_MM','/combination',num2str(c_comb)]);
        save(NMI_indi_path,'bestGamma','bestOmega','bestScore')
    else
        
        NMI_indi_path=fullfile(['/home/bianlb/Documents/CommuDetectLBM/HBI_CommuDetect/Results/real_Kong_Bayes_opt_MM','/combination',num2str(c_comb)]);
        save(NMI_indi_path,'bestGamma','bestOmega','bestScore')

    end
    close all

end
% -------------------------------------------------------------------------
% nested function

function val = objective_fun_group(x, A_global)

    gamma = x.gamma;
    omega = x.omega;
   
    N_subj = length(A_global);    

    score_sub_1 = zeros(N_subj,1);
    score_sub_2 = zeros(N_subj,1);
    
    z_1 = multilayer_modularity(A_global(:,1), gamma, omega);
    z_2 = multilayer_modularity(A_global(:,2), gamma, omega);
    z_3 = multilayer_modularity(A_global(:,3), gamma, omega);
    z_4 = multilayer_modularity(A_global(:,4), gamma, omega);
    
    % NMI
    for s = 1:N_subj  
        score_sub_1(s) = nmi(z_1(:,s),z_2(:,s));
        score_sub_2(s) = nmi(z_3(:,s),z_4(:,s));
    end
    val=-corr(score_sub_1,score_sub_2);
end







