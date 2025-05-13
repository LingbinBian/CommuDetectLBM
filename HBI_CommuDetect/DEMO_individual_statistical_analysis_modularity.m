% DEMO Statistical analysis: comparison between LBM and modularity
% individual-level analysis

% Version 1.0
% 7-April-2025
% Copyright (c) 2025, Lingbin Bian

clear
clc
close all

condition=3;
n_s=0.3162;
DIIV=30;
reso=[1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2];
reso_axis={'1','1.1','1.2','1.3','1.4','1.5','1.6','1.7','1.8','1.9','2'};
data1=zeros(100,length(reso));
data2=zeros(100,length(reso));

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



%age_range={'0-5 vs 6-11 Month', '3-8 vs 9-14 Month', '6-11 vs 12-17 Month','9-14 vs 15-23 Month','12-17 vs 18-29 Month','15-23 vs 24-36 Month','18-29 vs >36 Month'};

% data1=corre_hier;  % n by m matrix; n: number of samples, m: length of age range
% data2=corre_ave;   % n by m matrix; n: number of samples, m: length of age range



figure
% edge color
edgecolor1=[0,0,0]; % black color
edgecolor2=[0,0,0]; % black color

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

fillcolor1=[0,0,1]; % fillcolors = rand(24, 3);
fillcolor2=[1,0,0];

fillcolors=[repmat(fillcolor1,length(reso),1);repmat(fillcolor2,length(reso),1)];

position_1 = 0.8:1:(length(reso)-1+0.8);  % define position for first group boxplot
position_2 = 1.2:1:(length(reso)+0.2);  % define position for second group boxplot 
box_1 = boxplot(data1,'positions',position_1,'colors',edgecolor1,'width',0.2,'symbol','r+','outliersize',5);
hold on;
box_2 = boxplot(data2,'positions',position_2,'colors',edgecolor2,'width',0.2,'symbol','r+','outliersize',5);

boxobj = findobj(gca,'Tag','Box');
for j=1:length(boxobj)
    patch(get(boxobj(j),'XData'),get(boxobj(j),'YData'),fillcolors(j,:),'FaceAlpha',0.5,'LineWidth', 1);
end

set(gca,'XTick', 1:length(reso),'Xlim',[0 length(reso)+1]);
boxchi = get(gca, 'Children');

hLegend=legend([boxchi(1),boxchi(length(reso)+1)], ["Latent Block Model", "Modularity"],'location','southeast' );

ylim([0,1.5]); % range of y
set(gca,'xticklabel',reso_axis,'FontSize',12);
set(gca, 'linewidth', 1.2, 'fontsize', 16, 'fontname', 'times')   
xlabel('\gamma','fontsize',16)
ylabel('NMI','fontsize',16)
title(['Condition ',num2str(condition),', DIIV=',num2str(DIIV),])

set(gcf,'unit','centimeters','position',[6 10 18 12])
set(gca,'Position',[.22 .28 .75 .6]);

sig_x=zeros(length(reso),2);
sig_y=zeros(length(reso),1);

p_v=zeros(length(reso),1);
h=zeros(length(reso),1);
ci=cell(length(reso),1);
stats=cell(length(reso),1);

p=zeros(length(reso),1);
h1=zeros(length(reso),1);
ci1=cell(length(reso),1);

x=[0.8,1.2];

for i=1:length(reso)
    [h(i),p_v(i),ci{i},stats{i}] = vartest2(data1(:,i),data2(:,i));
    if p_v < 0.05
        [h1(i),p(i),ci1{i}] = ttest2(data1(:,i),data2(:,i),'Vartype','unequal');
    else
        [h1(i),p(i),ci1{i}] = ttest2(data1(:,i),data2(:,i));
    end
    sigline(x,max(max([data1(:,i),data2(:,i)])), p(i),p_v(i)); %
    x=x+1;
end

 
currentLabels = hLegend.String; % Get current labels (cell array)
Labels_remove={currentLabels{1},currentLabels{2}};
hLegend.String=Labels_remove;


data_path = fileparts(mfilename('fullpath'));

results_path=fullfile(data_path,['Results/synthetic_modularity/','DIIV',num2str(DIIV),'/n',num2str(n_s),'/statistical_analysis_individual_',num2str(condition)]);
save(results_path); 
saveas(gcf,['Results/synthetic_modularity/','DIIV',num2str(DIIV),'/n',num2str(n_s),'/statistical_analysis_individual_',num2str(condition),'_DIIV_',num2str(DIIV),'.fig'])
saveas(gcf,['Results/synthetic_modularity/','DIIV',num2str(DIIV),'/n',num2str(n_s),'/statistical_analysis_individual_',num2str(condition),'_DIIV_',num2str(DIIV),'.svg'])
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
    %print('not significance');
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
    fprintf('no significance ');
end


plot(x, [1;1]*y*1.1, '-k', 'LineWidth',1); % significance horizontal line
plot([1;1]*x(1), [y*1.05, y*1.1], '-k', 'LineWidth', 1); % significance vertical line
plot([1;1]*x(2), [y*1.05, y*1.1], '-k', 'LineWidth', 1); % significance vertical line

hold off
end




