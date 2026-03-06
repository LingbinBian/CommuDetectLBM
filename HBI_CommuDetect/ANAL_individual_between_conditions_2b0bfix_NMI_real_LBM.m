% MANIP (manipulate data) calculate the individual-level normalized mutual information (NMI) between conditions
%
% Version 1.1 
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

condition = {'2-back/fixation', '0-back/fixation'};

% Load task conditions
if atlas == 1
    load('Results/real_LBM/LR/grouplevel_data.mat');
else
    load('Results/real_Kong_LBM/LR/grouplevel_data.mat');
end
real_labels_L = z_g;

data1 = zeros(N_sub, 1); % For 2-back/fixation
data2 = zeros(N_sub, 1); % For 0-back/fixation

% Calculate NMI for the specified conditions
for i = 1:N_sub
    data1(i) = nmi(real_labels_L{i, 6}, real_labels_L{i, 4}); % 2-back vs fixation
    data2(i) = nmi(real_labels_L{i, 6}, real_labels_L{i, 5}); % 0-back vs fixation
end

% -------------------------------------------------------------------------
% Box plot (comparing two vectors of data)
figure
edgecolor1 = [0, 0, 0]; % black color
edgecolor2 = [0, 0, 0]; % black color

fillcolor1 = [0.5, 1, 0.82]; % Light green for 2-back/fixation
fillcolor2 = [0.94, 0.90, 0.55]; % Light yellow for 0-back/fixation
fillcolors = [fillcolor1; fillcolor2];

position_1 = 1;  
position_2 = 2;  

box_1 = boxplot(data1, 'positions', position_1, 'colors', edgecolor1, ...
                 'width', 0.5, 'symbol', 'r+', 'outliersize', 5);
hold on;
box_2 = boxplot(data2, 'positions', position_2, 'colors', edgecolor2, ...
                 'width', 0.5, 'symbol', 'r+', 'outliersize', 5);

boxobj = findobj(gca, 'Tag', 'Box');

for j = 1:length(boxobj)
    patch(get(boxobj(j), 'XData'), get(boxobj(j), 'YData'), fillcolors(j, :), ...
          'FaceAlpha', 0.5, 'LineWidth', 1);
end

set(gca, 'XTick', 1:length(condition), 'Xlim', [0 length(condition) + 1]);
ylim([0, 0.8]);
xlim([0.5, 2.5]);
set(gca, 'xticklabel', condition, 'FontSize', 12);
set(gca, 'linewidth', 1.2, 'fontsize', 16, 'fontname', 'times');   
xlabel('Condition', 'fontsize', 16);
ylabel('NMI', 'fontsize', 16);
title('Between Conditions');

set(gcf, 'unit', 'centimeters', 'position', [6 10 8 24]);
set(gca, 'Position', [.22 .28 .75 .6]);

% Perform t-test
[h, p_v, ci] = vartest2(data1, data2);
if p_v < 0.05
    [h_ttest, p_ttest, ci_ttest] = ttest2(data1, data2, 'Vartype', 'unequal');
else
    [h_ttest, p_ttest, ci_ttest] = ttest2(data1, data2);
end

% Signal significance line
x12 = [1, 2];
sigline(x12, max(max([data1, data2])), p_ttest, p_v); 

% Save results
data_path = fileparts(mfilename('fullpath'));
if atlas==1
    results_path=fullfile(data_path,'Results/real_LBM/LR/between_conditions_2b0bfix');
    save(results_path); 
    saveas(gcf,'Results/real_LBM/LR/between_conditions_2b0bfix.fig')
    saveas(gcf,'Results/real_LBM/LR/between_conditions_2b0bfix.svg')
else
    results_path=fullfile(data_path,'Results/real_LBM/LR/between_conditions_2b0bfix');
    save(results_path); 
    saveas(gcf,'Results/real_Kong_LBM/LR/between_conditions_2b0bfix.fig')
    saveas(gcf,'Results/real_Kong_LBM/LR/between_conditions_2b0bfix.svg')
end
% -------------------------------------------------------------------------
% sigline function
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