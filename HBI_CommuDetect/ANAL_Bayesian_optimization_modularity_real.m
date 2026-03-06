% ANAL Bayesian optimization for gamma in modularity model(real data)
%
% individual-level analysis
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
atlas=1; % SETTING 1: Schaefer 100, 2: Kong 200 

datatype=1;   % 1: real data, 0: synthetic data
subjects=load('subject.txt');
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
% Working memory fMRI data paradigm
% 1:2-back tool, 2:0-back body, 3:fixation, 4:2-back face, 5:0-back tool, 6:fixation, 7:2-back body,
% 8:2-back place, 9:fixation, 10:0-back face, 11:0-back place, 12:fixation

for c_comb=1:1
    switch c_comb
        case 1               % 2-back 
            c=[1,4,7,8];     % Tool & Face; Body & Place 
        case 2
            c=[1,7,4,8];     % Tool & Body; Face & Place 
        case 3
            c=[1,8,4,7];     % Tool & Place; Face & Body 
        case 4               % 0-back 
            c=[2,5,10,11];   % Tool & Face; Body & Place
        case 5
            c=[2,10,5,11];   % Tool & Body; Face & Place 
        case 6
            c=[2,11,5,10];   % Tool & Place; Face & Body 
    end
    fprintf('Combination: %d\n',c_comb)
  
    opt_parameter(group_adj,N_subj,atlas,c_comb,c);
   
end

% -------------------------------------------------------------------------
% nested function

function []=opt_parameter(group_adj,N_subj,atlas,c_comb,c)
    % Group-level Bayesian optimization (shared gamma)
    A_global = cell(N_subj,4);   % a combination of four working memory states (e.g. four 2-back states)

    for s = 1:N_subj
        A_global{s,1} = group_adj{s,1}{1,c(1)};
        A_global{s,2} = group_adj{s,1}{1,c(2)};
        A_global{s,3} = group_adj{s,1}{1,c(3)};
        A_global{s,4} = group_adj{s,1}{1,c(4)};
     
    end
    
    % Bayesian optimization    
    objective = @(x) objective_fun_group(x, A_global);
    
    gammaVar = optimizableVariable('gamma',[1.5 3.5]);
    
    results = bayesopt(objective, gammaVar,...
        'MaxObjectiveEvaluations',100,...
        'IsObjectiveDeterministic',true,...
        'AcquisitionFunctionName','expected-improvement-plus');
    
    bestGamma = results.XAtMinObjective.gamma;
    bestScore = results.MinObjective;
    
    fprintf('Best group-level gamma = %.4f\n', bestGamma);
    fprintf('Best objective value = %.6f\n', bestScore);
    
    
    % ---------------------------------------------------------------------
    % extract data
    
    gammaTrace = results.XTrace.gamma;
    scoreTrace = -results.ObjectiveTrace;
    
    bestGamma  = results.XAtMinObjective.gamma;
    bestScore  = max(scoreTrace);
    
    bestSoFar  = cummax(scoreTrace);
    
    % ---------------------------------------------------------------------
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

    if atlas==1
      
        para_path=fullfile(['/home/bianlb/Documents/CommuDetectLBM/HBI_CommuDetect/Results/real_Bayes_opt_M','/combination',num2str(c_comb)]);
        save(para_path,'bestGamma','bestScore')
    else
       
        para_path=fullfile(['/home/bianlb/Documents/CommuDetectLBM/HBI_CommuDetect/Results/real_Kong_Bayes_opt_M','/combination',num2str(c_comb)]);
        save(para_path,'bestGamma','bestScore')
    end

    close all

end
% -------------------------------------------------------------------------
% nested function

function val = objective_fun_group(x, A_global)
    % objective function
    gamma = x.gamma;

    N_subj = length(A_global);
    score_sub_1 = zeros(N_subj,1);
    score_sub_2 = zeros(N_subj,1);

    for s = 1:N_subj
        z_1 = modularity_und(A_global{s,1}, gamma);
        z_2 = modularity_und(A_global{s,2}, gamma);
        z_3 = modularity_und(A_global{s,3}, gamma);
        z_4 = modularity_und(A_global{s,4}, gamma);

        % normalized mutual information (NMI)
        score_sub_1(s) = nmi(z_1,z_2);
        score_sub_2(s) = nmi(z_3,z_4);
    end
    val=-corr(score_sub_1,score_sub_2);

end







