% MANIP (manipulate data) calculate the normalized mutual information (NMI) between estimation and ground truth
% group-level analysis using latent block model (LBM)
%
% Version 1.0 
% Copyright (c) 2025, Lingbin Bian
% 09-April-2025

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

 DIIV=30;


load(['Results/synthetic_LBM','/DIIV',num2str(DIIV),'/n',num2str(n_s),'/grouplevel_results.mat']);

% normalized mutual information
NMI=zeros(1,3);

NMI(1,1)=nmi(label_compare{1,1}(:,1),label_compare{1,1}(:,2));
NMI(1,2)=nmi(label_compare{1,2}(:,1),label_compare{1,2}(:,2));
NMI(1,3)=nmi(label_compare{1,3}(:,1),label_compare{1,3}(:,2));

data_path = fileparts(mfilename('fullpath'));
NMI_path=fullfile(data_path,['Results/synthetic_LBM','/DIIV',num2str(DIIV),'/n',num2str(n_s),'/NMI']);
save(NMI_path,'NMI')


