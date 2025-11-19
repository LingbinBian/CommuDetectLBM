% MANIP (manipulate data) calculate the between-subject normalized mutual information (NMI)

% Version 1.0 
% Copyright (c) 2025, Lingbin Bian
% 28-April-2025

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

% Load task condition data based on atlas selection
if atlas == 1
    load('Results/real_LBM/LR/grouplevel_data.mat');
else
    load('Results/real_Kong_LBM/LR/grouplevel_data.mat');
end
real_labels_L = z_g;

nmi_condition = zeros((N_sub * N_sub - N_sub) / 2, 11);

for c = 1:11
    labels_sub = zeros(N, N_sub);
    for i = 1:N_sub
        labels_sub(:, i) = real_labels_L{i, c};
    end
    nmi_condition(:, c) = betweensub_nmi(labels_sub);
end

length_nmi = (N_sub * N_sub - N_sub) / 2;
data = zeros(4 * length_nmi, 2);

% Collect data for boxplots
data(1:length_nmi, 1) = nmi_condition(:, 1);
data(length_nmi + 1:2 * length_nmi, 1) = nmi_condition(:, 4);
data(2 * length_nmi + 1:3 * length_nmi, 1) = nmi_condition(:, 7);
data(3 * length_nmi + 1:4 * length_nmi, 1) = nmi_condition(:, 8);

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
edgecolor = [0, 0, 0];
fillcolor1 = [0.78, 0.38, 0.08];
fillcolor2 = [0.53, 0.81, 0.92];
fillcolor3 = [1, 0.5, 0.5];
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
ylim([0, 0.8]); % range of y
xlim([0.5, 3.5])
set(gca, 'xticklabel', condition, 'FontSize', 12);
set(gca, 'linewidth', 1.2, 'fontsize', 16, 'fontname', 'times');
xlabel('Condition', 'fontsize', 16)
ylabel('NMI', 'fontsize', 16)
title('Between Subjects')

set(gcf, 'unit', 'centimeters', 'position', [6 10 8 24]);
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

if atlas==1
    results_path=fullfile(data_path,'Results/real_LBM/LR/between_subjects_2b0bfix');
    save(results_path); 
    saveas(gcf,'Results/real_LBM/LR/between_subjects_2b0bfix.fig')
    saveas(gcf,'Results/real_LBM/LR/between_subjects_2b0bfix.svg')
else
    results_path=fullfile(data_path,'Results/real_Kong_LBM/LR/between_subjects_2b0bfix');
    save(results_path); 
    saveas(gcf,'Results/real_Kong_LBM/LR/between_subjects_2b0bfix.fig')
    saveas(gcf,'Results/real_Kong_LBM/LR/between_subjects_2b0bfix.svg')
end
% ---

% -------------------------------------------------------------------------
% Significance line function
function sigline(x, y, p, p_v)
hold on
x = x';

if p < 0.001
    plot(mean(x), y * 1.15, '*k');
    plot(mean(x) - 0.16, y * 1.15, '*k');
    plot(mean(x) + 0.16, y * 1.15, '*k');
elseif (0.001 <= p) && (p < 0.01)
    plot(mean(x) - 0.08, y * 1.15, '*k');
    plot(mean(x) + 0.08, y * 1.15, '*k');
elseif (0.01 <= p) && (p < 0.05)
    plot(mean(x), y * 1.15, '*k');
end

if p_v < 0.001
    plot(mean(x), y * 1.2, 'pentagram', 'color', [0 0 0]);
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
function nmi_out = betweensub_nmi(labels_sub)
    [N_node, N_sub] = size(labels_sub);
    nmi_out = zeros((N_sub * N_sub - N_sub) / 2, 1);
    c = 1;
    for i = 1:N_sub
        for j = 1:N_sub
            if i < j
                nmi_out(c) = nmi(labels_sub(:, i), labels_sub(:, j));
                c = c + 1;
            end
        end
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