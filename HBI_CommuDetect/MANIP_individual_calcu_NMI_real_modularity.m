% MANIP (manipulate data) calculate the normalized mutual information (NMI)
% of the community memberships between 2-back tool and 2-back face
% individual-level analysis using modularity
%
% Version 1.0 
% Copyright (c) 2025, Lingbin Bian
% 21-April-2025

clear
clc
close all

atlas=2;
N_sub=100;
gamma=2;

% 1:2-back tool, 2:0-back body, 3:fixation, 4:2-back face, 5:0-back tool, 6:fixation, 7:2-back body,
% 8:2-back place, 9:fixation, 10:0-back face, 11:0-back place, 12:fixation

if atlas==1
    load(['Results/real_modularity/',num2str(gamma),'/LR/grouplevel_data.mat']);
else
    load(['Results/real_Kong_modularity/',num2str(gamma),'/LR/grouplevel_data.mat']);
end
real_labels_L=z_g;

%--------------------------------------------------------------------------
% Schaefer 100
% 2-back          % gamma=1.2    % gamma=1.3         % gamma=1.4        % gamma=1.5       % gamma=1.6        % gamma=1.7        % gamma=1.8         % gamma=1.9          % gamma=2
%c=[1,4,7,8];     % NaN          % r=0.2546          % r=0.2274         % r=0.2823        % r=0.2511         % r=0.2294         % r=0.2338          % r=0.2841           % r=0.2872
%c=[1,7,4,8];     % NaN          % r=0.2863          % r=0.2609         % r=0.2365        % r=0.1944         % r=0.1844         % r=0.1945          % r=0.2721           % r=0.2829
%c=[1,8,4,7];     % NaN          % r=0.3061          % r=0.3254         % r=0.3285        % r=0.3443         % r=0.3530         % r=0.3913          % r=0.3899           % r=0.4269
                  % NaN          % r_mean=0.2823     % r_mean=0.2712    % r_mean=0.2824   % r_mean=0.2633    % r_mean=0.2556    % r_mean=0.2732     % r_mean=0.3154      % r_mean=0.3323
%--------------------------------------------------------------------------
% 0-back          % gamma=1.2    % gamma=1.3         % gamma=1.4        % gamma=1.5       % gamma=1.6        % gamma=1.7        % gamma=1.8         % gamma=1.9           % gamma=2           
%c=[2,5,10,11];   % NaN          % r=0.2248          % r=0.1956         % r=0.1755        % r=0.2001         % r=0.2664         % r=0.2624          % r=0.2791            % r=0.2438     
%c=[2,10,5,11];   % NaN          % r=0.4380          % r=0.4153         % r=0.4506        % r=0.4085         % r=0.4348         % r=0.4042          % r=0.3407            % r=0.3019           
%c=[2,11,5,10];   % NaN          % r=0.3513          % r=0.3483         % r=0.3172        % r=0.3245         % r=0.3240         % r=0.2681          % r=0.3246            % r=0.3079       
                  % NaN          % r_mean=0.3380     % r_mean=0.3197    % r_mean=0.3144   % r_mean=0.3110    % r_mean=0.3417    % r_mean=0.3116     % r_mean=0.3148       % r_mean=0.2845

%--------------------------------------------------------------------------
% Kong 200 
% 2-back         % gamma=1      % gamma=1.5       % gamma=2
%c=[1,4,7,8];    % NaN          % r=0.3438        % r=0.3141
%c=[1,7,4,8];    % NaN          % r=0.2462        % r=0.2638
%c=[1,8,4,7];    % NaN          % r=0.3232        % r=0.3472
                                % r_mean=0.3044   % r_mean=0.3084
%--------------------------------------------------------------------------
% 0-back         % gamma=1     % gamma=1.5        % gamma=2
%c=[2,5,10,11];  % NaN         % r=0.2708         % r=0.2592 
%c=[2,10,5,11];  % NaN         % r=0.4918         % r=0.4087
%c=[2,11,5,10];  % NaN         % r=0.3690         % r=0.3534
                               % r_mean=0.3772    % r_mean=0.3404
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

r=corr(half_1,half_2)
% 
% data_path = fileparts(mfilename('fullpath'));
% NMI_indi_path=fullfile(data_path,['Results/real_modularity/',num2str(gamma),'/NMI_indi_real_modularity']);
% save(NMI_indi_path,'NMI_indi_real_modularity')


