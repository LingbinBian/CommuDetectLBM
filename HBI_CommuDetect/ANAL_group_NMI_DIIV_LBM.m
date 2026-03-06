% MANIP NMI with different levels of DIIV at the group level
%
% Version 1.0 
% Copyright (c) 2025, Lingbin Bian
% 31-March-2025

clear
clc
close all

% signal to noise ratio SNR=10log10(sigma_signal^2/sigma_noise^2);
% sigma_signal=1, sigma_noise=n_s
 n_s=0.3162;  % 10dB
% n_s=0.5623;  % 5dB
% n_s=1;    % 0dB  
% n_s=1.7783;  % -5dB
% n_s=3.1623;  % -10dB
% n_s=5.6234;  % -15dB
% n_s=10;      % -20dB
% n_s=17.7828; % -25dB
DIIV=[10,20,30];

% n_s=[10, 5.6234, 3.1623, 1.7783, 1, 0.5623, 0.3162];

% snr_vector=[-20, -15, -10, -5, 0, 5, 10];

NMI_vs_DIIV=zeros(3,length(DIIV));
for i=1:length(DIIV)
    load(['Results/synthetic_LBM','/DIIV',num2str(DIIV(i)),'/n',num2str(n_s),'/NMI.mat'])
    NMI_vs_DIIV(:,i)=NMI;
end

%colorvector=[1,1,0;0.78,0.38,0.08;0,0,1;1,0,0;0,1,0;0,0.5,0;0.5,0.5,0;1,0.5,0.5];
colorvector=[1,1,0;0.78,0.38,0.08;0,0,1;1,0,0;0,1,0;0,0.5,0;0.5,0.5,0;1,0.5,0.5];


for i=1:3
    figure
    plot(DIIV(:),NMI_vs_DIIV(i,:),'--ks',...
    'LineWidth',1.2,...
    'MarkerSize',7,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',colorvector(2,:));
    hold on  
    
    title(['Condition ',num2str(i),', DIIV=',num2str(DIIV(i)),])
    set(gca,'box','on')
    set(gca,'fontsize',16)
    xlabel('SNR(dB)','fontsize',16)
    ylabel('NMI','fontsize',16)
    xticks([10,20,30]);
    % set(gcf,'unit','centimeters','position',[6 10 8 5])
    % set(gca,'Position',[.25 .3 .55 .55]);
   
    set(gcf,'unit','centimeters','position',[6 10 10 10])
    set(gca,'Position',[.22 .28 .75 .6]);

    %set(gca,'xtick',0:2:20)
    set(gca, 'linewidth', 1.2, 'fontsize', 16, 'fontname', 'times')
    saveas(gcf,['Results/synthetic_LBM/DIIV',num2str(DIIV(i)),'/Group_LBM_SNR_10dB_',num2str(i),'_DIIV_',num2str(DIIV(i)),'.fig'])
    saveas(gcf,['Results/synthetic_LBM/DIIV',num2str(DIIV(i)),'/Group_LBM_SNR_10dB_',num2str(i),'_DIIV_',num2str(DIIV(i)),'.svg'])
end




