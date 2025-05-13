% MANIP (manipulate data) calculate the normalized mutual information (NMI)
% of the community memberships between 2-back tool and 2-back face
% individual-level analysis using latent block model (LBM)
%
% Version 1.0 
% Copyright (c) 2025, Lingbin Bian
% 20-April-2025

clear
clc
close all

N_sub=100;
atlas=2;
% Task conditions

% 1:2-back tool, 2:0-back body, 3:fixation, 4:2-back place, 5:0-back tool, 6:fixation, 7:2-back body,
% 8:2-back place, 9:fixation, 10:0-back face, 11:0-back place, 12:fixation

if atlas==1
    load('Results/real_LBM/LR/grouplevel_data.mat');
else
    load('Results/real_Kong_LBM/LR/grouplevel_data.mat');
end
real_labels_L=z_g;

%--------------------------------------------------------------------------
% Schaefer 100
% 2-back
%c=[1,4,7,8];  % r=0.5166
%c=[1,7,4,8];  % r=0.4915
%c=[1,8,4,7];  % r=0.4707
               % r_mean=0.4929
%--------------------------------------------------------------------------
% 0-back
%c=[2,5,10,11];  % r=0.3233
%c=[2,10,5,11];  % r=0.3995
%c=[2,11,5,10];  % r=0.4662
                % r_mean=0.3963


%--------------------------------------------------------------------------
% Kong 200
% 2-back
%c=[1,4,7,8];  % r=0.5147
%c=[1,7,4,8];  % r=0.5282
%c=[1,8,4,7];  % r=0.4462
               % r_mean=0.4964
%--------------------------------------------------------------------------
% 0-back
%c=[2,5,10,11];  % r=0.4411
%c=[2,10,5,11];  % r=0.3981
%c=[2,11,5,10];  % r=0.3868
                 % r_mean=0.4087
% normalized mutual information

half_1=zeros(N_sub,1);
half_2=zeros(N_sub,1);


for i=1:N_sub
        label_temporal=zeros(N,2);
        label_temporal(:,1)=real_labels_L{i,c(1)};
        label_temporal(:,2)=real_labels_L{i,c(2)};
        label_temporal=labelswitch(label_temporal);
        real_labels_L{i,c(1)}=label_temporal(:,1);
        real_labels_L{i,c(2)}=label_temporal(:,2);
        half_1(i,1)=nmi(real_labels_L{i,c(1)},real_labels_L{i,c(2)});
        label_temporal=zeros(N,2);
end
   
for i=1:N_sub
        label_temporal=zeros(N,2);
        label_temporal(:,1)=real_labels_L{i,c(3)};
        label_temporal(:,2)=real_labels_L{i,c(4)};
        label_temporal=labelswitch(label_temporal);
        real_labels_L{i,c(3)}=label_temporal(:,1);
        real_labels_L{i,c(4)}=label_temporal(:,2);
        half_2(i,1)=nmi(real_labels_L{i,c(3)},real_labels_L{i,c(4)});
        label_temporal=zeros(N,2);
end

r_correlation=corr(half_1,half_2)

% data_path = fileparts(mfilename('fullpath'));
% NMI_indi_path=fullfile(data_path,'Results/real/NMI_indi_real_LBM');
% save(NMI_indi_path,'NMI_indi_real_LBM')



