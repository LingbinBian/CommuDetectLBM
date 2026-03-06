% ANAL Bayesian optimization for LBM (real data)
%
%
% Version 1.0
% 25-Feb-2026
% Copyright (c) 2026, Lingbin Bian
% -------------------------------------------------------------------------
clear
clc
close all
% -------------------------------------------------------------------------
% paralell computing
ncore = str2double(getenv('SLURM_CPUS_PER_TASK'));
if isempty(gcp('nocreate'))
    parpool('local', ncore);
end

% Load data
% Data type
atlas=1; % SETTING 1: Schaefer 100, 2: Kong 200 
fprintf('Atlas: %d\n',atlas)
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

for c_comb=1:6
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
    % Group-level Bayesian optimization
    A_global = cell(N_subj,4);   % a combination of four working memory states (e.g. four 2-back states)

    for s = 1:N_subj       
        A_global{s,1} = group_adj{s,1}{1,c(1)};
        A_global{s,2} = group_adj{s,1}{1,c(2)};
        A_global{s,3} = group_adj{s,1}{1,c(3)};
        A_global{s,4} = group_adj{s,1}{1,c(4)};
     
    end
    
    % Bayesian optimization    
    objective = @(x) objective_fun_group(x, A_global);
    
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

    if atlas==1
       
        para_path=fullfile(['/home/bianlb/Documents/CommuDetectLBM/HBI_CommuDetect/Results/real_Bayes_opt_LBM','/combination',num2str(c_comb)]);
        save(para_path,'bestNu','bestRho','bestKappa_sq','bestK_max','bestScore')
    else
        
        para_path=fullfile(['/home/bianlb/Documents/CommuDetectLBM/HBI_CommuDetect/Results/real_Kong_Bayes_opt_LBM','/combination',num2str(c_comb)]);
        save(para_path,'bestNu','bestRho','bestKappa_sq','bestK_max','bestScore')
    end

    close all

end
% -------------------------------------------------------------------------
% nested function

function val = objective_fun_group(x, A_global)

    nu       = x.nu;
    rho      = x.rho;
    kappa_sq = x.kappa_sq;
    K_max    = x.K_max;

    N_subj = length(A_global);
    score_sub_1 = zeros(N_subj,1);
    score_sub_2 = zeros(N_subj,1);

    parfor s = 1:N_subj
        
        A1 = A_global{s,1};
        A2 = A_global{s,2};
        A3 = A_global{s,3};
        A4 = A_global{s,4};
        fprintf('Subject: %d\n',s)
        z_1 = CommuDetectLBM(A1,K_max,nu,rho,kappa_sq);
        z_2 = CommuDetectLBM(A2,K_max,nu,rho,kappa_sq);
        z_3 = CommuDetectLBM(A3,K_max,nu,rho,kappa_sq);
        z_4 = CommuDetectLBM(A4,K_max,nu,rho,kappa_sq);

        score_sub_1(s) = nmi(z_1,z_2);
        score_sub_2(s) = nmi(z_3,z_4);
    end

    val = -corr(score_sub_1,score_sub_2);

end







