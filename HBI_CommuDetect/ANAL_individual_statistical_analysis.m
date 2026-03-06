% DEMO Statistical analysis: comparison between LBM, modularity, and multilayer modularity
% individual-level analysis (added multiple comparison correction: FDR)
%
% Version 1.1 
% 17-Nov-2025
% Copyright (c) 2025, Lingbin Bian

clear
clc
close all

condition=3; % for each data segment
n_s=0.3162;
DIIV=30;
reso=[1,1.5,2]; % modularity resolution gamma
reso_axis={'1','1.5','2'};
data1=zeros(100,length(reso));
data2=zeros(100,length(reso));
data3=zeros(100,length(reso));
data4=zeros(100,length(reso));
data5=zeros(100,length(reso));

load(['Results/synthetic_LBM/','DIIV',num2str(DIIV),'/n',num2str(n_s),'/NMI_indi_LBM.mat']);

for i=1:length(reso)
    data1(:,i)=NMI_indi_LBM(:,condition);
end


for i=1:length(reso)
    load(['Results/synthetic_modularity/','DIIV',num2str(DIIV),'/n',num2str(n_s),'/',num2str(reso(i)),'/NMI_indi_modularity.mat']);
    NMI_indi_temporal=NMI_indi_modularity(:,condition);
    data2(:,i)=NMI_indi_temporal;
    NMI_indi_temporal=zeros(100,1);
end

omega=0.01;
for i=1:length(reso)
    load(['Results/synthetic_multilayer_modularity_coupling',num2str(omega),'/','DIIV',num2str(DIIV),'/n',num2str(n_s),'/',num2str(reso(i)),'/NMI_indi_multilayer_modularity.mat']);
    NMI_indi_temporal=NMI_indi_multilayer_modularity(:,condition);
    data3(:,i)=NMI_indi_temporal;
    NMI_indi_temporal=zeros(100,1);
end
omega=0.1;
for i=1:length(reso)
    load(['Results/synthetic_multilayer_modularity_coupling',num2str(omega),'/','DIIV',num2str(DIIV),'/n',num2str(n_s),'/',num2str(reso(i)),'/NMI_indi_multilayer_modularity.mat']);
    NMI_indi_temporal=NMI_indi_multilayer_modularity(:,condition);
    data4(:,i)=NMI_indi_temporal;
    NMI_indi_temporal=zeros(100,1);
end
omega=0.5;
for i=1:length(reso)
    load(['Results/synthetic_multilayer_modularity_coupling',num2str(omega),'/','DIIV',num2str(DIIV),'/n',num2str(n_s),'/',num2str(reso(i)),'/NMI_indi_multilayer_modularity.mat']);
    NMI_indi_temporal=NMI_indi_multilayer_modularity(:,condition);
    data5(:,i)=NMI_indi_temporal;
    NMI_indi_temporal=zeros(100,1);
end

% figure and plotting
figure
% edge color
edgecolor1=[0,0,0]; % black color
edgecolor2=[0,0,0]; % black color
edgecolor3=[0,0,0]; % black color
edgecolor4=[0,0,0]; % black color
edgecolor5=[0,0,0]; % black color

fillcolor5=[0,0,1]; % fillcolors = rand(24, 3);
fillcolor4=[1,0,0];
fillcolor3=[0.75 0.75 0.75];
fillcolor2=[0.78,0.38,0.08;];
fillcolor1=[0.5,0.5,0];

fillcolors=[repmat(fillcolor1,length(reso),1);repmat(fillcolor2,length(reso),1);repmat(fillcolor3,length(reso),1);repmat(fillcolor4,length(reso),1);repmat(fillcolor5,length(reso),1)];

position_1 = 0.7:1:(length(reso)-1+0.7);  % define position for first group boxplot
position_2 = 0.85:1:(length(reso)-1+0.85);  % define position for first group boxplot
position_3 = 1.0:1:(length(reso));  % define position for second group boxplot 
position_4 = 1.15:1:(length(reso)+0.15);  % define position for second group boxplot 
position_5 = 1.3:1:(length(reso)+0.3);  % define position for second group boxplot 

box_1 = boxplot(data1,'positions',position_1,'colors',edgecolor1,'width',0.075,'symbol','r+','outliersize',5);
hold on;
box_2 = boxplot(data2,'positions',position_2,'colors',edgecolor2,'width',0.075,'symbol','r+','outliersize',5);
hold on
box_3 = boxplot(data3,'positions',position_3,'colors',edgecolor3,'width',0.075,'symbol','r+','outliersize',5);
hold on
box_4 = boxplot(data4,'positions',position_4,'colors',edgecolor4,'width',0.075,'symbol','r+','outliersize',5);
hold on
box_5 = boxplot(data5,'positions',position_5,'colors',edgecolor5,'width',0.075,'symbol','r+','outliersize',5);
hold on

boxobj = findobj(gca,'Tag','Box');
for j=1:length(boxobj)
    patch(get(boxobj(j),'XData'),get(boxobj(j),'YData'),fillcolors(j,:),'FaceAlpha',0.5,'LineWidth', 1);
end

set(gca,'XTick', 1:length(reso),'Xlim',[0 length(reso)+1]);
boxchi = get(gca, 'Children');

hLegend=legend([boxchi(1),boxchi(4),boxchi(7),boxchi(10),boxchi(13)], ["Latent Block Model","Modularity","Multilayer Modularity \omega=0.01","Multilayer Modularity \omega=0.1","Multilayer Modularity \omega=0.5"],'location','southeast','FontSize', 10 );
xlim([0.5,3.5]); % range of y

ylim([0.3,1.5]); % range of y
set(gca,'xticklabel',reso_axis,'FontSize',12);
set(gca, 'linewidth', 1.2, 'fontsize', 16, 'fontname', 'times')   
xlabel('\gamma','fontsize',16)
ylabel('NMI','fontsize',16)
title(['Experiment ',num2str(condition),', DIIV=',num2str(DIIV),])

set(gcf,'unit','centimeters','position',[6 10 20 16])
set(gca,'Position',[.22 .28 .75 .6]);

% -------------------------------------------------------------------------
% New statistics block: collect p-values, run FDR (BH), and plot corrected siglines
% -------------------------------------------------------------------------

% Preallocate containers
num_reso = length(reso);
num_comp = 4 * num_reso; % 4 comparisons per gamma
p_all = zeros(1,num_comp);      % t-test p-values
p_v_all = zeros(1,num_comp);    % vartest p-values (variance test)
positions_all = cell(1,num_comp);
y_all = zeros(1,num_comp);
comp_names = cell(1,num_comp);

idx = 0;

% 1) LBM vs Modularity
x = [0.7,0.85];
for i=1:num_reso
    idx = idx + 1;
    [h_v,p_v] = vartest2(data1(:,i),data2(:,i));
    if p_v < 0.05
        [h_t,p] = ttest2(data1(:,i),data2(:,i),'Vartype','unequal');
    else
        [h_t,p] = ttest2(data1(:,i),data2(:,i));
    end
    p_all(idx) = p;
    p_v_all(idx) = p_v;
    positions_all{idx} = x;
    y_all(idx) = max(max([data1(:,i),data2(:,i)]));
    comp_names{idx} = sprintf('LBM vs Modularity (gamma=%g)', reso(i));
    x = x + 1;
end

% 2) LBM vs Multilayer (omega = 0.01)
x = [0.7,1.0];
for i=1:num_reso
    idx = idx + 1;
    [h_v,p_v] = vartest2(data1(:,i),data3(:,i));
    if p_v < 0.05
        [h_t,p] = ttest2(data1(:,i),data3(:,i),'Vartype','unequal');
    else
        [h_t,p] = ttest2(data1(:,i),data3(:,i));
    end
    p_all(idx) = p;
    p_v_all(idx) = p_v;
    positions_all{idx} = x;
    y_all(idx) = 1.12 * max(max([data1(:,i),data3(:,i)]));
    comp_names{idx} = sprintf('LBM vs MMult w=0.01 (gamma=%g)', reso(i));
    x = x + 1;
end

% 3) LBM vs Multilayer (omega = 0.1)
x = [0.7,1.15];
for i=1:num_reso
    idx = idx + 1;
    [h_v,p_v] = vartest2(data1(:,i),data4(:,i));
    if p_v < 0.05
        [h_t,p] = ttest2(data1(:,i),data4(:,i),'Vartype','unequal');
    else
        [h_t,p] = ttest2(data1(:,i),data4(:,i));
    end
    p_all(idx) = p;
    p_v_all(idx) = p_v;
    positions_all{idx} = x;
    y_all(idx) = 1.25 * max(max([data1(:,i),data4(:,i)]));
    comp_names{idx} = sprintf('LBM vs MMult w=0.1 (gamma=%g)', reso(i));
    x = x + 1;
end

% 4) LBM vs Multilayer (omega = 0.5)
x = [0.7,1.3];
for i=1:num_reso
    idx = idx + 1;
    [h_v,p_v] = vartest2(data1(:,i),data5(:,i));
    if p_v < 0.05
        [h_t,p] = ttest2(data1(:,i),data5(:,i),'Vartype','unequal');
    else
        [h_t,p] = ttest2(data1(:,i),data5(:,i));
    end
    p_all(idx) = p;
    p_v_all(idx) = p_v;
    positions_all{idx} = x;
    y_all(idx) = 1.35 * max(max([data1(:,i),data5(:,i)]));
    comp_names{idx} = sprintf('LBM vs MMult w=0.5 (gamma=%g)', reso(i));
    x = x + 1;
end

% Perform FDR correction on t-test p-values and vartest p-values
if exist('mafdr','file') == 2
    % Use MATLAB's mafdr if available
    p_fdr = mafdr(p_all, 'BHFDR', true);
    p_v_fdr = mafdr(p_v_all, 'BHFDR', true);
else
    % Fallback to Benjamini-Hochberg implementation (returns corrected p-values)
    p_fdr = benjamini_hochberg(p_all);
    p_v_fdr = benjamini_hochberg(p_v_all);
end

% Draw significance markers using corrected p-values
for k = 1:length(p_all)
    % Use corrected p-values for plotting
    sigline(positions_all{k}, y_all(k), p_fdr(k), p_v_fdr(k));
end

% Print a summary table
fprintf('\nComparison results (original p, FDR-corrected p):\n');
for k = 1:length(p_all)
    fprintf('%2d) %s\n   t-test p = %.4g, FDR p = %.4g; vartest p = %.4g, FDR p = %.4g\n', ...
        k, comp_names{k}, p_all(k), p_fdr(k), p_v_all(k), p_v_fdr(k));
end

% Update legend labels (unchanged behavior)
currentLabels = hLegend.String; % Get current labels (cell array)
Labels_remove={currentLabels{1},currentLabels{2},currentLabels{3},currentLabels{4},currentLabels{5}};
hLegend.String=Labels_remove;

data_path = fileparts(mfilename('fullpath'));

% optionally save results
results_path=fullfile(data_path,['Results/synthetic_LBM/','DIIV',num2str(DIIV),'/n',num2str(n_s),'/statistical_analysis_individual_LBM_M_MM',num2str(condition)]);
save(results_path); 
saveas(gcf,['Results/synthetic_LBM/','DIIV',num2str(DIIV),'/n',num2str(n_s),'/statistical_analysis_individual_LBM_M_MM',num2str(condition),'_DIIV_',num2str(DIIV),'.fig'])
saveas(gcf,['Results/synthetic_LBM/','DIIV',num2str(DIIV),'/n',num2str(n_s),'/statistical_analysis_individual_LBM_M_MM',num2str(condition),'_DIIV_',num2str(DIIV),'.svg'])

% -------------------------------------------------------------------------
% sigline (unchanged, uses supplied p and p_v values)
function sigline(x, y, p, p_v)
hold on
x = x';

if p<0.001
    plot(mean(x),       y*1.075, '*k')          % the sig star sign
    plot(mean(x)- 0.06, y*1.075, '*k')          % the sig star sign
    plot(mean(x)+ 0.06, y*1.075, '*k')          % the sig star sign

elseif (0.001<=p)&&(p<0.01)
    plot(mean(x)- 0.03, y*1.075, '*k')         % the sig star sign
    plot(mean(x)+ 0.03, y*1.075, '*k')          % the sig star sign

elseif (0.01<=p)&&(p<0.05)
    plot(mean(x), y*1.075, '*k')               % the sig star sign
else
    % not significant for t-test (after correction)
end

if p_v<0.001
    plot(mean(x),       y*1.105, 'pentagram','color',[0 0 0])          % the sig star sign
    plot(mean(x)- 0.06, y*1.105, 'pentagram','color',[0 0 0])          % the sig star sign
    plot(mean(x)+ 0.06, y*1.105, 'pentagram','color',[0 0 0])          % the sig star sign

elseif (0.001<=p_v)&&(p_v<0.01)
    plot(mean(x)- 0.03, y*1.105, 'pentagram','color',[0 0 0])         % the sig star sign
    plot(mean(x)+ 0.03, y*1.105, 'pentagram','color',[0 0 0])         % the sig star sign

elseif (0.01<=p_v)&&(p_v<0.05)
    plot(mean(x), y*1.105, 'pentagram','color',[0 0 0])               % the sig star sign
else
    % not significant for vartest (after correction)
end

plot(x, [1;1]*y*1.05, '-k', 'LineWidth',0.5); % significance horizontal line
plot([1;1]*x(1), [y*1.025, y*1.05], '-k', 'LineWidth', 0.5); % significance vertical line
plot([1;1]*x(2), [y*1.025, y*1.05], '-k', 'LineWidth', 0.5); % significance vertical line

hold off
end

% -------------------------------------------------------------------------
% Benjamini-Hochberg FDR implementation (returns adjusted p-values)
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
            p_adj_sorted(i) = min(p_sorted(i) * m / i, p_adj_sorted(i+1));
        end
    end
    
    % Cap at 1
    p_adj_sorted = min(p_adj_sorted, 1);
    
    % Unsort to original order
    p_adj = zeros(size(p));
    p_adj(sort_idx) = p_adj_sorted;
    p_adj = reshape(p_adj, size(p));
end

