% DEMO Statistical analysis: NMI(LBM) vs SNR
% individual-level analysis

% Version 1.0
% 18-April-2025
% Copyright (c) 2025, Lingbin Bian

clear
clc
close all

condition=1;

% signal to noise ratio SNR=10log10(sigma_signal^2/sigma_noise^2);
% sigma_signal=1, sigma_noise:n_s
% n_s=0.3162;  % 10dB
% n_s=0.5623;  % 5dB
% n_s=1;       % 0dB
% n_s=1.7783;  % -5dB
% n_s=3.1623;  % -10dB
% n_s=5.6234;  % -15dB
% n_s=10;      % -20dB
% n_s=17.7828; % -25dB

DIIV=10;
n_s=[10,5.6234,3.1623,1.7783,1,0.5623,0.3162]; % standard deviation of noise
snr={'-20','-15','-10','-5','0','5','10'}; % SNR in dB
data=zeros(100,length(n_s)); 

for i=1:length(n_s)
    load(['Results/synthetic_LBM/','DIIV',num2str(DIIV),'/n',num2str(n_s(i)),'/NMI_indi_LBM.mat']);
    NMI_indi_temporal=NMI_indi_LBM(:,condition);
    data(:,i)=NMI_indi_temporal;
    NMI_indi_temporal=zeros(100,1);
end

figure
% edge color

edgecolor=[0,0,0]; % black color

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

fillcolor1=[0.5,0.54,0.53]; % fillcolors = rand(24, 3);
fillcolor2=[0.5,0.54,0.53];

fillcolors=[repmat(fillcolor1,length(n_s),1);repmat(fillcolor2,length(n_s),1)];

position = 1:1:(length(n_s));  % define position for second group boxplot 
%box_1 = boxplot(data1,'positions',position_1,'colors',edgecolor1,'width',0.2,'symbol','r+','outliersize',5);
%hold on;
box = boxplot(data,'positions',position,'colors',edgecolor,'width',0.2,'symbol','r+','outliersize',5);

boxobj = findobj(gca,'Tag','Box');
for j=1:length(boxobj)
    patch(get(boxobj(j),'XData'),get(boxobj(j),'YData'),fillcolors(j,:),'FaceAlpha',0.5,'LineWidth', 1);
end

set(gca,'XTick', 1:length(n_s),'Xlim',[0 length(n_s)+1]);
boxchi = get(gca, 'Children');

%hLegend=legend([boxchi(1),boxchi(length(n_s)+1)], ["Latent Block Model", "Modularity"],'location','southeast' );

ylim([0,1]); % range of y
set(gca,'xticklabel',snr,'FontSize',12);
set(gca, 'linewidth', 1.2, 'fontsize', 16, 'fontname', 'times')   
xlabel('SNR (dB)','fontsize',16)
ylabel('NMI','fontsize',16)
title(['Condition ',num2str(condition),', DIIV=',num2str(DIIV),])

set(gcf,'unit','centimeters','position',[6 10 10 10])
set(gca,'Position',[.22 .28 .75 .6]);



mean_data=mean(data);
variance_data=var(data);

data_path = fileparts(mfilename('fullpath'));

results_path=fullfile(data_path,['Results/synthetic_LBM/','DIIV',num2str(DIIV),'/statistical_analysis_LBM_SNR_',num2str(condition)]);
save(results_path); 
saveas(gcf,['Results/synthetic_LBM/','DIIV',num2str(DIIV),'/statistical_analysis_LBM_SNR_',num2str(condition),'_DIIV_',num2str(DIIV),'.fig'])
saveas(gcf,['Results/synthetic_LBM/','DIIV',num2str(DIIV),'/statistical_analysis_LBM_SNR_',num2str(condition),'_DIIV_',num2str(DIIV),'.svg'])
% -------------------------------------------------------------------------
% sigline

