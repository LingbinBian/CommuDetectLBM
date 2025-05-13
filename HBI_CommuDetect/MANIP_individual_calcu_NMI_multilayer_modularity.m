% MANIP (manipulate data) calculate the normalized mutual information (NMI)
% between estimation and ground truth
% individual-level analysis using multilayer modularity
%
% Version 1.0 
% Copyright (c) 2025, Lingbin Bian
% 12-April-2025

clear
clc
close all

% signal to noise ratio SNR=10log10 (sigma_signal^2/sigma_noise^2);
% sigma_signal=1, sigma_noise=n_s
 n_s=0.3162;  % 10dB
% n_s=0.5623;  % 5dB
% n_s=1;    % 0dB  
% n_s=1.7783;  % -5dB
% n_s=3.1623;  % -10dB
% n_s=5.6234;  % -15dB
% n_s=10;      % -20dB
% n_s=17.7828; % -25dB

% modularity resolution gamma=1:0.1:2

gamma=2;
omega=0.5;

DIIV=30;
N_sub=100;
N_state=3;

load(['Results/synthetic_multilayer_modularity_coupling',num2str(omega),'/','DIIV',num2str(DIIV),'/n',num2str(n_s),'/',num2str(gamma),'/grouplevel_data.mat']);

% normalized mutual information

NMI_indi_multilayer_modularity=zeros(100,N_state);

for i=1:N_sub
    for j=1:N_state
        NMI_indi_multilayer_modularity(i,j)=nmi(z_g{i,j},true_latent_sub{i,1}(:,j));
    end
end

data_path = fileparts(mfilename('fullpath'));

NMI_indi_path=fullfile(data_path,['Results/synthetic_multilayer_modularity_coupling',num2str(omega),'/','DIIV',num2str(DIIV),'/n',num2str(n_s),'/',num2str(gamma),'/NMI_indi_multilayer_modularity']);
save(NMI_indi_path,'NMI_indi_multilayer_modularity')


