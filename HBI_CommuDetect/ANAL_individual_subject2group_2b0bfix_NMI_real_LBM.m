% MANIP (manipulate data) calculate the subject-specific to group normalized mutual information (NMI)
% 2-back, 0-back, and fixation

% Version 1.0 
% Copyright (c) 2025, Lingbin Bian
% 30-April-2025

clear
clc
close all

N_sub = 100;
atlas = 2;
if atlas == 1
    N_node = 100;
else
    N_node = 200;
end

condition = {'2-back', '0-back', 'fixation'};

% Task conditions
if atlas == 1
    load('Results/real_LBM/LR/grouplevel_data.mat');
else
    load('Results/real_Kong_LBM/LR/grouplevel_data.mat');
end
real_labels_L = z_g;

if atlas == 1
    load('Results/real_LBM/LR/grouplevel_results.mat');
else
    load('Results/real_Kong_LBM/LR/grouplevel_results.mat');
end

labels_group = label_group_esti;
nmi_condition = zeros(N_sub, 11);

for c = 1:11
    labels_sub = zeros(N, N_sub);
    
    for i = 1:N_sub
        labels_sub(:, i) = real_labels_L{i, c};
    end
    nmi_condition(:, c) = subject2group(labels_sub, labels_group(:, c));
    labels_sub = zeros(N, N_sub);
end

length_nmi = N_sub;
data = zeros(4 * length_nmi, 2);

% 2-back
data(1:length_nmi, 1) = nmi_condition(:, 1);
data(length_nmi + 1:2 * length_nmi, 1) = nmi_condition(:, 4);
data(2 * length_nmi + 1:3 * length_nmi, 1) = nmi_condition(:, 7);
data(3 * length_nmi + 1:4 * length_nmi, 1) = nmi_condition(:, 8);

% 0-back
data(1:length_nmi, 2) = nmi_condition(:, 2);
data(length_nmi + 1:2 * length_nmi, 2) = nmi_condition(:, 5);
data(2 * length_nmi + 1:3 * length_nmi, 2) = nmi_condition(:, 10);
data(3 * length_nmi + 1:4 * length_nmi, 2) = nmi_condition(:, 11);

data1 = data(:, 1);
data2 = data(:, 2);

data3 = zeros(3 * length_nmi, 1);
data3(1:length_nmi, 1) = nmi_condition(:, 3);
data3(length_nmi + 1:2 * length_nmi, 1) = nmi_condition(:, 6);
data3(2 * length_nmi + 1:3 * length_nmi, 1) = nmi_condition(:, 9);

% -------------------------------------------------------------------------
% Box plot (comparing two vectors of data)
figure
edgecolor1 = [0, 0, 0];
edgecolor2 = [0, 0, 0];
edgecolor3 = [0, 0, 0];

fillcolor1 = [0.78, 0.38, 0.08]; 
fillcolor2 = [0.53, 0.81, 0.92]; 
fillcolor3 = [1, 0.5, 0.5]; 
fillcolors = [fillcolor1; fillcolor2; fillcolor3];

position_1 = 1;  
position_2 = 2;  
position_3 = 3;  

box_1 = boxplot(data1, 'positions', position_1, 'colors', edgecolor1, 'width', 0.5, 'symbol', 'r+', 'outliersize', 5);
hold on;
box_2 = boxplot(data2, 'positions', position_2, 'colors', edgecolor2, 'width', 0.5, 'symbol', 'r+', 'outliersize', 5);
hold on;
box_3 = boxplot(data3, 'positions', position_3, 'colors', edgecolor3, 'width', 0.5, 'symbol', 'r+', 'outliersize', 5);

boxobj = findobj(gca, 'Tag', 'Box');
for j = 1:length(boxobj)
    patch(get(boxobj(j), 'XData'), get(boxobj(j), 'YData'), fillcolors(j, :), 'FaceAlpha', 0.5, 'LineWidth', 1);
end

set(gca, 'XTick', 1:length(condition), 'Xlim', [0 length(condition) + 1]);
ylim([0, 0.8]);
xlim([0.5, 3.5]);
set(gca, 'xticklabel', condition, 'FontSize', 12);
set(gca, 'linewidth', 1.2, 'fontsize', 16, 'fontname', 'times');   
xlabel('Condition', 'fontsize', 16);
ylabel('NMI', 'fontsize', 16);
title('Subject-specific to Group');

set(gcf, 'unit', 'centimeters', 'position', [6 10 8 24]);
set(gca, 'Position', [.22 .28 .7 .6]);

% Perform t-tests
p_values = zeros(3, 1);
x12 = [1, 2];
[h_12, p_v_12] = vartest2(data1, data2);
if p_v_12 < 0.05
    [h1_12, p_12, ci1_12] = ttest2(data1, data2, 'Vartype', 'unequal');
else
    [h1_12, p_12, ci1_12] = ttest2(data1, data2);
end
p_values(1) = p_12;

x13 = [1, 3];
[h_13, p_v_13] = vartest2(data1, data3);
if p_v_13 < 0.05
    [h1_13, p_13, ci1_13] = ttest2(data1, data3, 'Vartype', 'unequal');
else
    [h1_13, p_13, ci1_13] = ttest2(data1, data3);
end
p_values(2) = p_13;

x23 = [2, 3];
[h_23, p_v_23] = vartest2(data2, data3);
if p_v_23 < 0.05
    [h1_23, p_23, ci1_23] = ttest2(data2, data3, 'Vartype', 'unequal');
else
    [h1_23, p_23, ci1_23] = ttest2(data2, data3);
end
p_values(3) = p_23;

% Apply Benjamini-Hochberg FDR control
adjusted_p_values = benjamini_hochberg(p_values);

% Plot significance lines based on corrected p-values
if adjusted_p_values(1) < 0.05
    sigline(x12, max(max([data1, data2])), p_12, p_v_12); 
end

if adjusted_p_values(2) < 0.05
    sigline(x13, 0.1 + max([max(data1), max(data3)]), p_13, p_v_13); 
end

if adjusted_p_values(3) < 0.05
    sigline(x23, 0.2 + max([max(data2), max(data3)]), p_23, p_v_23); 
end

% Save results
data_path = fileparts(mfilename('fullpath'));
if atlas == 1
    results_path = fullfile(data_path, 'Results/real_LBM/LR/subject2group_2b0bfix');
    save(results_path); 
    saveas(gcf, 'Results/real_LBM/LR/subject2group_2b0bfix.fig');
    saveas(gcf, 'Results/real_LBM/LR/subject2group_2b0bfix.svg');
else
    results_path = fullfile(data_path, 'Results/real_Kong_LBM/LR/subject2group_2b0bfix');
    save(results_path); 
    saveas(gcf, 'Results/real_Kong_LBM/LR/subject2group_2b0bfix.fig');
    saveas(gcf, 'Results/real_Kong_LBM/LR/subject2group_2b0bfix.svg');
end

% -------------------------------------------------------------------------
% sigline function
function sigline(x, y, p, p_v)
hold on
x = x';

if p < 0.001
    plot(mean(x),       y * 1.15, '*k');         
    plot(mean(x) - 0.16, y * 1.15, '*k');         
    plot(mean(x) + 0.16, y * 1.15, '*k');         
elseif (0.001 <= p) && (p < 0.01)
    plot(mean(x) - 0.08, y * 1.15, '*k');         
    plot(mean(x) + 0.08, y * 1.15, '*k');         
elseif (0.01 <= p) && (p < 0.05)
    plot(mean(x), y * 1.15, '*k');               
end

if p_v < 0.001
    plot(mean(x),       y * 1.2, 'pentagram', 'color', [0 0 0]);         
    plot(mean(x) - 0.16, y * 1.2, 'pentagram', 'color', [0 0 0]);         
    plot(mean(x) + 0.16, y * 1.2, 'pentagram', 'color', [0 0 0]);         
elseif (0.001 <= p_v) && (p_v < 0.01)
    plot(mean(x) - 0.08, y * 1.2, 'pentagram', 'color', [0 0 0]);         
    plot(mean(x) + 0.08, y * 1.2, 'pentagram', 'color', [0 0 0]);         
elseif (0.01 <= p_v) && (p_v < 0.05)
    plot(mean(x), y * 1.2, 'pentagram', 'color', [0 0 0]);               
end

plot(x, [1; 1] * y * 1.12, '-k', 'LineWidth', 1); 
plot([1; 1] * x(1), [y * 1.1, y * 1.12], '-k', 'LineWidth', 1); 
plot([1; 1] * x(2), [y * 1.1, y * 1.12], '-k', 'LineWidth', 1); 

hold off
end

% -------------------------------------------------------------------------
% Nested function for NMI calculation
function nmi_out = subject2group(labels_sub, labels_group)

[N_node, N_sub] = size(labels_sub);
nmi_out = zeros(N_sub, 1);

for i = 1:N_sub
    nmi_out(i) = nmi(labels_sub(:, i), labels_group);
end
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