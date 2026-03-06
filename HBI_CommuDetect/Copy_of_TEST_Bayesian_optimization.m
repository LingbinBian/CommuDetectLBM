% Test Bayesian optimization

% simulate functional connectivity
clear; clc; rng(1);

N = 120;                 % node
K = 4;                   % true number of communities
T = 300;                 % time course

% true community labels
true_label = repelem(1:K, N/K)';

% generate block covariance
Sigma = 0.2*ones(N);
for k = 1:K
    idx = true_label==k;
    Sigma(idx,idx) = 0.6;
end

Sigma = Sigma + 0.4*eye(N);

% 生成模拟 fMRI 时间序列
X = mvnrnd(zeros(N,1), Sigma, T);

% 计算相关矩阵
A = corr(X);

% 保证非负（modularity_und要求）
A(A<0) = 0;


% -------------------------------------------------------------------------
% Bayesian optimization

A_global = A;   % global variable (adjacency matrix)

objective = @(x) objective_fun(x, A_global);

gammaVar = optimizableVariable('gamma',[0.5 2.5]);

results = bayesopt(objective, gammaVar,...
    'MaxObjectiveEvaluations',30,...
    'IsObjectiveDeterministic',false,...
    'AcquisitionFunctionName','expected-improvement-plus');

bestGamma = results.XAtMinObjective.gamma

% -------------------------------------------------------------------------

figure;
plot(results.XTrace.gamma, -results.ObjectiveTrace,'o-');
xlabel('\gamma');
ylabel('Stability');
title('Bayesian optimization trace');

figure;
plot(-results.ObjectiveTrace,'o-','LineWidth',1.5);
xlabel('Iteration');
ylabel('Stability');
title('Optimization Progress');
grid on;


figure;
scatter(results.XTrace.gamma, -results.ObjectiveTrace,60,'filled');
xlabel('\gamma');
ylabel('Stability');
title('Explored \gamma vs Stability');
grid on;


bestSoFar = cummax(-results.ObjectiveTrace);

figure;
plot(bestSoFar,'LineWidth',1.5);
xlabel('Iteration');
ylabel('Best Stability So Far');
title('Best Objective Value');
grid on;


% -------------------------------------------------------------------------
% nested function
function val = objective_fun(x, A)

    gamma = x.gamma;

    stability = modularity_stability(gamma, A);

    val = -stability;   % bayesopt (default minimum value)
end

function score = modularity_stability(gamma, A)

    nRuns = 20;   % repetition time
    labels = [];

    for r = 1:nRuns
        z = modularity_und(A, gamma);
        labels(:,r) = z;
    end

    % 计算 pairwise NMI
    cnt = 1;
    vals = [];

    for i = 1:nRuns
        for j = i+1:nRuns
            vals(cnt) = nmi(labels(:,i), labels(:,j));
            cnt = cnt + 1;
        end
    end

    score = mean(vals);  % stability
end

