% MANIP (manipulate data) calculate the individual-level normalized mutual information (NMI) between conditions
%
% Version 1.0 
% Copyright (c) 2025, Lingbin Bian
% 30-April-2025

clear
clc
close all

N_sub=100;
atlas=2;
if atlas==1
    N_node=100;
else
    N_node=200;
end

condition={'2-back/0-back','2-back/fixation','0-back/fixation'};

% Task conditions

% 1:2-back tool, 2:0-back body, 3:fixation, 4:2-back place, 5:0-back tool, 6:fixation, 7:2-back body,
% 8:2-back place, 9:fixation, 10:0-back face, 11:0-back place, 12:fixation

if atlas==1
    load('Results/real_LBM/LR/grouplevel_data.mat');
else
    load('Results/real_Kong_LBM/LR/grouplevel_data.mat');
end
real_labels_L=z_g;


data1=zeros(N_sub,1);
data2=zeros(N_sub,1);
data3=zeros(N_sub,1);
% for i=1:N_sub
%     data1(i)=nmi(real_labels_L{i,1},real_labels_L{i,2});
%     data2(i)=nmi(real_labels_L{i,1},real_labels_L{i,3});
%     data3(i)=nmi(real_labels_L{i,2},real_labels_L{i,3});
% end

% 4:2-back place, 5:0-back tool, 6:fixation
for i=1:N_sub
    data1(i)=nmi(real_labels_L{i,4},real_labels_L{i,5});
    data2(i)=nmi(real_labels_L{i,4},real_labels_L{i,6});
    data3(i)=nmi(real_labels_L{i,5},real_labels_L{i,6});
end
% 
% for i=1:N_sub
%     data1(i)=nmi(real_labels_L{i,8},real_labels_L{i,10});
%     data2(i)=nmi(real_labels_L{i,8},real_labels_L{i,11});
%     data3(i)=nmi(real_labels_L{i,10},real_labels_L{i,11});
% end





% -------------------------------------------------------------------------
% blox plot (comparing two vectors of data)
figure
% edge color
edgecolor1=[0,0,0]; % black color
edgecolor2=[0,0,0]; % black color
edgecolor3=[0,0,0]; % black color

%fillcolor

% 0.75 0.75 0.75
% 1,0.87,0.68
% 1,1,0;
% 0.78,0.38,0.08;
% 0,0,1;
% 1,0,0;
% 0,1,0;
% 0,0.5,0;
% 0.5,0.5,0;
% 1,0.5,0.5;

fillcolor1=[0.87,0.63,0.87]; % 
fillcolor2=[0.5,1,0.82]; % 
fillcolor3=[0.94,0.90,0.55];
fillcolors=[fillcolor1;fillcolor2;fillcolor3];

position_1 = 1;  % define position for first group boxplot
position_2 = 2;  % define position for second group boxplot 
position_3 = 3;  % define position for second group boxplot 

box_1 = boxplot(data1,'positions',position_1,'colors',edgecolor1,'width',0.5,'symbol','r+','outliersize',5);
hold on;
box_2 = boxplot(data2,'positions',position_2,'colors',edgecolor2,'width',0.5,'symbol','r+','outliersize',5);
hold on;
box_3 = boxplot(data3,'positions',position_3,'colors',edgecolor3,'width',0.5,'symbol','r+','outliersize',5);

boxobj = findobj(gca,'Tag','Box');

for j=1:length(boxobj)
    patch(get(boxobj(j),'XData'),get(boxobj(j),'YData'),fillcolors(j,:),'FaceAlpha',0.5,'LineWidth', 1);
end

set(gca,'XTick', 1:length(condition),'Xlim',[0 length(condition)+1]);
boxchi = get(gca, 'Children');

%hLegend=legend([boxchi(1),boxchi(length(condition)+1)], ["2-back", "0-back"],'location','southeast' );

ylim([0,0.9]); % range of y
xlim([0.5,3.5])
set(gca,'xticklabel',condition,'FontSize',12);
set(gca, 'linewidth', 1.2, 'fontsize', 16, 'fontname', 'times')   
xlabel('Condition','fontsize',16)
ylabel('NMI','fontsize',16)
title('Between Conditions')

set(gcf,'unit','centimeters','position',[6 10 8 24])
set(gca,'Position',[.22 .28 .75 .6]);

x12=[1,2];

[h_12,p_v_12,ci_12,stats_12] = vartest2(data1,data2);
if p_v_12 < 0.05
    [h1_12,p_12,ci1_12] = ttest2(data1,data2,'Vartype','unequal');
else
    [h1_12,p_12,ci1_12] = ttest2(data1,data2);
end
sigline(x12,max(max([data1,data2])),p_12,p_v_12); %


x13=[1,3];

[h_13,p_v_13,ci_13,stats_13] = vartest2(data1,data3);
if p_v_13 < 0.05
    [h1_13,p_13,ci1_13] = ttest2(data1,data3,'Vartype','unequal');
else
    [h1_13,p_13,ci1_13] = ttest2(data1,data3);
end
sigline(x13,0.1+max([max(data1),max(data3)]),p_13,p_v_13); %



x23=[2,3];

[h_23,p_v_23,ci_23,stats_23] = vartest2(data2,data3);
if p_v_23 < 0.05
    [h1_23,p_23,ci1_23] = ttest2(data2,data3,'Vartype','unequal');
else
    [h1_23,p_23,ci1_23] = ttest2(data2,data3);
end
sigline(x23,0.21+max([max(data2),max(data3)]),p_23,p_v_23); %


% currentLabels = hLegend.String; % Get current labels (cell array)
% Labels_remove={currentLabels{1},currentLabels{2}};
% hLegend.String=Labels_remove;



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
% sigline
function sigline(x, y, p, p_v)
hold on
x = x';

if p<0.001
    plot(mean(x),       y*1.15, '*k')          % the sig star sign
    plot(mean(x)- 0.16, y*1.15, '*k')          % the sig star sign
    plot(mean(x)+ 0.16, y*1.15, '*k')          % the sig star sign

elseif (0.001<=p)&&(p<0.01)
    plot(mean(x)- 0.08, y*1.15, '*k')         % the sig star sign
    plot(mean(x)+ 0.08, y*1.15, '*k')         % the sig star sign

elseif (0.01<=p)&&(p<0.05)
    plot(mean(x), y*1.15, '*k')               % the sig star sign
else
    %fprintf('no significance');
end

if p_v<0.001
    plot(mean(x),       y*1.2, 'pentagram','color',[0 0 0])          % the sig star sign
    plot(mean(x)- 0.16, y*1.2, 'pentagram','color',[0 0 0])          % the sig star sign
    plot(mean(x)+ 0.16, y*1.2, 'pentagram','color',[0 0 0])          % the sig star sign

elseif (0.001<=p_v)&&(p_v<0.01)
    plot(mean(x)- 0.08, y*1.2, 'pentagram','color',[0 0 0])         % the sig star sign
    plot(mean(x)+ 0.08, y*1.2, 'pentagram','color',[0 0 0])         % the sig star sign

elseif (0.01<=p_v)&&(p_v<0.05)
    plot(mean(x), y*1.2, 'pentagram','color',[0 0 0])               % the sig star sign
else
    %fprintf('no significance');
end


plot(x, [1;1]*y*1.1, '-k', 'LineWidth',1); % significance horizontal line
plot([1;1]*x(1), [y*1.05, y*1.1], '-k', 'LineWidth', 1); % significance vertical line
plot([1;1]*x(2), [y*1.05, y*1.1], '-k', 'LineWidth', 1); % significance vertical line

hold off
end



% -------------------------------------------------------------------------
% nested function



