% MANIP NMI with different levels of SNR at the group level
%
% Version 1.0 
% Copyright (c) 2025, Lingbin Bian
% 31-March-2025

clear
clc
close all

% signal to noise ratio SNR=10log10(sigma_signal^2/sigma_noise^2);
% sigma_signal=1, sigma_noise=n_s
% n_s=0.3162;  % 10dB
% n_s=0.5623;  % 5dB
% n_s=1;    % 0dB  
 n_s=1.7783;  % -5dB
% n_s=3.1623;  % -10dB
% n_s=5.6234;  % -15dB
% n_s=10;      % -20dB
% n_s=17.7828; % -25dB
DIIV=10;

n_s=[10, 5.6234, 3.1623, 1.7783, 1, 0.5623, 0.3162];

snr_vector=[-20, -15, -10, -5, 0, 5, 10];

NMI_vs_snr=zeros(3,length(n_s));
NMI_consensus_vs_snr=zeros(3,length(n_s));
NMI_majorityvote_vs_snr=zeros(3,length(n_s));
for i=1:length(n_s)
    load(['Results/synthetic_LBM','/DIIV',num2str(DIIV),'/n',num2str(n_s(i)),'/NMI.mat'])
    NMI_vs_snr(:,i)=NMI;
end

% for i=1:length(n_s)
%     load(['Results/synthetic_LBM','/DIIV',num2str(DIIV),'/n',num2str(n_s(i)),'/NMI_consensus.mat'])
%     NMI_consensus_vs_snr(:,i)=NMI;
% end

for i=1:length(n_s)
    load(['Results/synthetic_LBM','/DIIV',num2str(DIIV),'/n',num2str(n_s(i)),'/NMI_majorityvote.mat'])
    NMI_majorityvote_vs_snr(:,i)=NMI;
end

%colorvector=[1,1,0;0.78,0.38,0.08;0,0,1;1,0,0;0,1,0;0,0.5,0;0.5,0.5,0;1,0.5,0.5];
colorvector=[1,1,0;0.78,0.38,0.08;0,0,1;1,0,0;0,1,0;0,0.5,0;0.5,0.5,0;1,0.5,0.5];


for i=1:3
    figure 
    plot(snr_vector(:),NMI_majorityvote_vs_snr(i,:),'--ks',...
    'LineWidth',1.8,...
    'MarkerSize',8,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',colorvector(1,:));
    hold on
    % plot(snr_vector(:),NMI_consensus_vs_snr(i,:),'--ks',...
    % 'LineWidth',1.8,...
    % 'MarkerSize',8,...
    % 'MarkerEdgeColor','k',...
    % 'MarkerFaceColor',colorvector(3,:));
    % hold on  
    plot(snr_vector(:),NMI_vs_snr(i,:),'--ks',...
    'LineWidth',1.8,...
    'MarkerSize',8,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',colorvector(4,:));
    hold on 
    
    title(['Condition ',num2str(i),', DIIV=',num2str(DIIV),])
    set(gca,'box','on')
    set(gca,'fontsize',16)
    xlabel('SNR(dB)','fontsize',16)
    ylabel('NMI','fontsize',16)
    xticks([-20, -15, -10, -5, 0, 5, 10]);
    % set(gcf,'unit','centimeters','position',[6 10 8 5])
    % set(gca,'Position',[.25 .3 .55 .55]);
    legend('MV','CD','location','southeast')
    set(gcf,'unit','centimeters','position',[6 10 10 10])
    set(gca,'Position',[.22 .28 .75 .6]);

    %set(gca,'xtick',0:2:20)
    set(gca, 'linewidth', 1.2, 'fontsize', 16, 'fontname', 'times')
    saveas(gcf,['Results/synthetic_LBM/DIIV',num2str(DIIV),'/Group_LBM_SNR_',num2str(i),'_DIIV_',num2str(DIIV),'.fig'])
    saveas(gcf,['Results/synthetic_LBM/DIIV',num2str(DIIV),'/Group_LBM_SNR_',num2str(i),'_DIIV_',num2str(DIIV),'.svg'])
end




