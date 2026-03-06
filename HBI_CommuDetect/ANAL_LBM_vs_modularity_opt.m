% LBM vs (Multilayer) modularity (synthetic data)
%
% Version 1.0 
% Copyright (c) 2026, Lingbin Bian
% 22-Feb-2026

clear
clc
close all

N_sub = 100;
atlas = 1;

DIIV =30;
switch DIIV
    case 10
        n_diiv=1;
    case 20
        n_diiv=2;
    otherwise
        n_diiv=3;
end
state=3;
n_s=0.3162;

if atlas == 1
    N_node = 100;
else
    N_node = 200;
end

condition = {'LBM', 'M', 'MM'};

% Load task conditions
if atlas == 1
    %load(['Results/synthetic_LBM/','DIIV',num2str(DIIV),'/n',num2str(n_s),'/NMI_indi_LBM.mat']);

    load('Results/synthetic_Bayes_opt_LBM/infer_opt.mat');
    nmi_score_LBM=nmi_score;
   
    load('Results/synthetic_Bayes_opt_M/infer_opt.mat');
    nmi_score_M=nmi_score;

    load('Results/synthetic_Bayes_opt_MM/infer_opt.mat');
    nmi_score_MM=nmi_score;


else
    load('Results/real_Kong_LBM/LR/grouplevel_data.mat');
end


data1 = zeros(N_sub, 1); % LBM
data2 = zeros(N_sub, 1); % modularity
data3 = zeros(N_sub, 1); % multilayer modularity

% Calculate NMI for the specified method
for i = 1:N_sub
    data1(i) = nmi_score_LBM{state,n_diiv}(i,1); 
    data2(i) = nmi_score_M{state,n_diiv}(i,1); 
    data3(i) = nmi_score_MM{state,n_diiv}(i,1); 

end

% -------------------------------------------------------------------------
% Box plot (comparing two vectors of data)
figure
edgecolor = [0, 0, 0];
% fillcolor1 = [0.5,0.5,0];
% fillcolor2 = [0.53, 0.81, 0.92];
% fillcolor3 = [1, 0.5, 0.5];

fillcolor3=[0,0,1]; % fillcolors = rand(24, 3);
fillcolor2=[1,0,0];
fillcolor1=[0.75 0.75 0.75];
% fillcolor2=[0.78,0.38,0.08;];
% fillcolor1=[0.5,0.5,0];

fillcolors = [fillcolor1; fillcolor2; fillcolor3];

% Define positions for boxplots
box_1 = boxplot(data1, 'positions', 1, 'colors', edgecolor, 'width', 0.5, 'symbol', 'r+', 'outliersize', 5);
hold on;
box_2 = boxplot(data2, 'positions', 2, 'colors', edgecolor, 'width', 0.5, 'symbol', 'r+', 'outliersize', 5);
box_3 = boxplot(data3, 'positions', 3, 'colors', edgecolor, 'width', 0.5, 'symbol', 'r+', 'outliersize', 5);

boxobj = findobj(gca, 'Tag', 'Box');
for j = 1:length(boxobj)
    patch(get(boxobj(j), 'XData'), get(boxobj(j), 'YData'), fillcolors(j,:), 'FaceAlpha', 0.5, 'LineWidth', 1);
end

set(gca, 'XTick', 1:length(condition), 'Xlim', [0 length(condition) + 1]);
ylim([0.5, 1.5]); % range of y
xlim([0.5, 3.5])
set(gca, 'xticklabel', condition, 'FontSize', 12);
set(gca, 'linewidth', 1.2, 'fontsize', 16, 'fontname', 'times');
xlabel('', 'fontsize', 16)
ylabel('NMI', 'fontsize', 16)
title(['Experiment ',num2str(state),', DIIV ',num2str(DIIV)])

set(gcf, 'unit', 'centimeters', 'position', [6 10 8 16]);
set(gca, 'Position', [.22 .28 .75 .6]);

% Perform t-tests
p_values = zeros(3, 1);
[~, p_values(1)] = ttest2(data1, data2); % 2-back vs 0-back
[~, p_values(2)] = ttest2(data1, data3); % 2-back vs fixation
[~, p_values(3)] = ttest2(data2, data3); % 0-back vs fixation

% Apply Benjamini-Hochberg FDR control
adjusted_p_values = benjamini_hochberg(p_values);

% Significance determination after correction
significant_flags = adjusted_p_values < 0.05;

% Plot significance lines based on corrected p-values
if significant_flags(1)
    sigline([1, 2], max(max([data1, data2])) + 0.05, p_values(1), p_values(1));
end

if significant_flags(2)
    sigline([1, 3], 0.1 + max([max(data1), max(data3)]), p_values(2), p_values(2));
end

if significant_flags(3)
    sigline([2, 3], 0.2 + max([max(data2), max(data3)]), p_values(3), p_values(3));
end

data_path = fileparts(mfilename('fullpath'));


saveas(gcf,['Results/LBM_vs_modularity','_experiment',num2str(state),'_DIIV',num2str(DIIV),'.fig'])
saveas(gcf,['Results/LBM_vs_modularity','_experiment',num2str(state),'_DIIV',num2str(DIIV),'.svg'])



% -------------------------------------------------------------------------
% Significance line function
function sigline(x, y, p, p_v)
    hold on
    x = x';
    
    mark_size=5;
    if p < 0.001
        plot(mean(x), y * 1.1, '*k','MarkerSize', mark_size);
        plot(mean(x) - 0.16, y * 1.1, '*k','MarkerSize', mark_size);
        plot(mean(x) + 0.16, y * 1.1, '*k','MarkerSize', mark_size);
    elseif (0.001 <= p) && (p < 0.01)
        plot(mean(x) - 0.08, y * 1.1, '*k','MarkerSize', mark_size);
        plot(mean(x) + 0.08, y * 1.1, '*k','MarkerSize', mark_size);
    elseif (0.01 <= p) && (p < 0.05)
        plot(mean(x), y * 1.1, '*k','MarkerSize', mark_size);
    end
    
    if p_v < 0.001
        plot(mean(x), y * 1.12, 'pentagram', 'color', [0 0 0],'MarkerSize', mark_size);
        plot(mean(x) - 0.16, y * 1.12, 'pentagram', 'color', [0 0 0],'MarkerSize', mark_size);
        plot(mean(x) + 0.16, y * 1.12, 'pentagram', 'color', [0 0 0],'MarkerSize', mark_size);
    elseif (0.001 <= p_v) && (p_v < 0.01)
        plot(mean(x) - 0.08, y * 1.12, 'pentagram', 'color', [0 0 0],'MarkerSize', mark_size);
        plot(mean(x) + 0.08, y * 1.12, 'pentagram', 'color', [0 0 0],'MarkerSize', mark_size);
    elseif (0.01 <= p_v) && (p_v < 0.05)
        plot(mean(x), y * 1.12, 'pentagram', 'color', [0 0 0],'MarkerSize', mark_size);
    end
    
    plot(x, [1; 1] * y * 1.085, '-k', 'LineWidth', 1);
    plot([1; 1] * x(1), [y * 1.06, y * 1.085], '-k', 'LineWidth', 1);
    plot([1; 1] * x(2), [y * 1.06, y * 1.085], '-k', 'LineWidth', 1);
    
    hold off
end



% -------------------------------------------------------------------------
% Benjamini-Hochberg FDR implementation function
function p_adj = benjamini_hochberg(p)
    % Input: p: vector of p-values
    % Output: p_adj: BH-adjusted p-values (same length), monotonic non-decreasing
    if isempty(p)
        p_adj = p;
        return;
    end
    
    % Sort p-values
    [p_sorted, sort_idx] = sort(p(:));
    m = length(p_sorted);
    
    % Allocate adjusted p-value array
    p_adj_sorted = zeros(size(p_sorted));
    
    % Compute BH adjusted p-values (step-up procedure)
    for i = m:-1:1
        if i == m
            p_adj_sorted(i) = p_sorted(i);
        else
            p_adj_sorted(i) = min(p_sorted(i) * m / i, p_adj_sorted(i + 1));
        end
    end
    
    % Cap at 1
    p_adj_sorted = min(p_adj_sorted, 1);
    
    % Unsort to original order
    p_adj = zeros(size(p));
    p_adj(sort_idx) = p_adj_sorted;
    p_adj = reshape(p_adj, size(p));
end