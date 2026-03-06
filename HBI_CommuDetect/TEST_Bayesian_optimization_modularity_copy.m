% Bayesian optimization for modularity
%
% individual-level analysis
%
% Version 1.0
% 6-March-2025
% Copyright (c) 2025, Lingbin Bian
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
    atlas=2;
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
   % n_s=0.3162;  % 10dB
   % n_s=0.5623;  % 5dB
    n_s=1;      
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

M=10; % margin size
if atlas==1
    N=100; 
else
    N=200; 
end


% -------------------------------------------------------------------------
% Bayesian optimization

A_global = group_adj{8,1}{1,1};   % global variable (adjacency matrix)
true_label=true_latent_sub{8,1}(:,1);
K=numel(unique(true_label));

% -------------------------------------------------------------------------
% Pre-scan gamma to estimate adaptive lambda

gamma_grid = linspace(1,2.5, 50);

NMIvals = zeros(size(gamma_grid));
Kvals   = zeros(size(gamma_grid));

for i = 1:length(gamma_grid)
    z = modularity_und(A_global, gamma_grid(i));
    NMIvals(i) = nmi(z, true_label);
    Kvals(i)   = numel(unique(z));
end

range_NMI = max(NMIvals) - min(NMIvals);
range_K   = max(abs(Kvals - K));

lambda = 0.5 * range_NMI / (range_K + eps);

fprintf('Adaptive lambda = %.4f\n', lambda);

% -------------------------------------------------------------------------


objective = @(x) objective_fun(x, A_global,true_label, lambda);

gammaVar = optimizableVariable('gamma',[1 2.5]);

results = bayesopt(objective, gammaVar,...
    'MaxObjectiveEvaluations',30,...
    'IsObjectiveDeterministic',true,...
    'AcquisitionFunctionName','expected-improvement-plus');

bestGamma = results.XAtMinObjective.gamma;
fprintf('Best gamma = %.4f\n', bestGamma);

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


function val = objective_fun(x, A, true_label, lambda)

    gamma = x.gamma;

    z = modularity_und(A, gamma);

    % NMI
    nmi_score = nmi(z, true_label);

    % Community number penalty
    K_gamma = numel(unique(z));
    K_true  = numel(unique(true_label));

    penalty = abs(K_gamma - K_true);

    score = nmi_score - lambda * penalty;

    val = -score;  % bayesopt minimizes
end





