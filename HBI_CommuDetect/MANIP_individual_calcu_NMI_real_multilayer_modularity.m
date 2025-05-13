% MANIP (manipulate data) calculate the normalized mutual information (NMI)
% of the community memberships between 2-back tool and 2-back face
% individual-level analysis using multilayer modularity
%
% Version 1.0 
% Copyright (c) 2025, Lingbin Bian
% 24-April-2025

clear
clc
close all

N_sub=100;
gamma=2; %(1:1.5)
omega=0.01;
atlas=2;
% 1:2-back tool, 2:0-back body, 3:fixation, 4:2-back face, 5:0-back tool, 6:fixation, 7:2-back body,
% 8:2-back place, 9:fixation, 10:0-back face, 11:0-back place, 12:fixation
if atlas==1
    load(['Results/real_multilayer_modularity_coupling',num2str(omega),'/',num2str(gamma),'/LR/grouplevel_data.mat']);
else
    load(['Results/real_Kong_multilayer_modularity_coupling',num2str(omega),'/',num2str(gamma),'/LR/grouplevel_data.mat']);
end
real_labels_L=z_g;

%--------------------------------------------------------------------------
% Schaefer 100
% omega=0.01
% 2-back         % gamma=1      % gamma=1.1       % gamma=1.2           % gamma=1.3          % gamma=1.4         % gamma=1.5         % gamma=1.6        % gamma=1.7         % gamma=1.8         % gamma=1.9        % gamma=2      
%c=[1,4,7,8];    % r=*          % r=*             % r=0.3026            % r=0.2907           % r=0.2945          % r=0.2511          % r=0.2306         % r=0.2427          % r=0.2504          % r=0.2304         % r=0.2442
%c=[1,7,4,8];    % r=*          % r=*             % r=0.3318            % r=0.2425           % r=0.2304          % r=0.2040          % r=0.2060         % r=0.2661          % r=0.2633          % r=0.2510         % r=0.2505 
%c=[1,8,4,7];    % r=*          % r=*             % r=0.3097            % r=0.2333           % r=0.3326          % r=0.3071          % r=0.3602         % r=0.3808          % r=0.3929          % r=0.3804         % r=0.3763
                 % r_mean=*     % r_mean=*        % r_mean=0.3147       % r_mean=0.2555      % r_mean=0.2858     % r_mean=0.2541     % r_mean=0.2656    % r_mean=0.2965     % r_mean=0.3022     % r_mean=0.2873    % r_mean=0.2903 
%--------------------------------------------------------------------------
% Schaefer 100
% omega=0.01
% 0-back         % gamma=1      % gamma=1.1       % gamma=1.2           % gamma=1.3          % gamma=1.4         % gamma=1.5         % gamma=1.6       % gamma=1.7          % gamma=1.8         % gamma=1.9        % gamma=2
%c=[2,5,10,11];  % r=*          % r=*             % r=0.3231            % r=0.2224           % r=0.1732          % r=0.2037          % r=0.2440        % r=0.2154           % r=0.2656          % r=0.2355         % r=0.2840
%c=[2,10,5,11];  % r=*          % r=*             % r=0.4620            % r=0.4330           % r=0.4193          % r=0.3871          % r=0.3624        % r=0.3799           % r=0.3673          % r=0.3416         % r=0.3560
%c=[2,11,5,10];  % r=*          % r=*             % r=0.3785            % r=0.3153           % r=0.3058          % r=0.2715          % r=0.2851         % r=0.2469          % r=0.2972          % r=0.2658         % r=0.2906
                 % r_mean=*     % r_mean=*        % r_mean=0.3879       % r_mean=0.3236      % r_mean=0.2994     % r_mean=0.2874     % r_mean=0.2972   % r_mean=0.2807      % r_mean=0.3100     % r_mean=0.2810    % r_mean=0.3102

%--------------------------------------------------------------------------
% Schaefer 100
% omega=0.1
% 2-back         % gamma=1      % gamma=1.1       % gamma=1.2        % gamma=1.3     % gamma=1.4       % gamma=1.5       % gamma=1.6     %         
%c=[1,4,7,8];    % r=*          % r=#             % r=#              % r=#           % r=#             % r=#             % r=#           %    
%c=[1,7,4,8];    % r=*          % r=#             % r=#              % r=#           % r=#             % r=#             % r=#           % 
%c=[1,8,4,7];    % r=*          % r=#             % r=#              % r=#           % r=#             % r=#             % r=#           %       
                 % r_mean=*     % r_mean=#        % r_mean=#         % r_mean=#      % r_mean=#        % r_mean=#        % r_mean=#      %  
%--------------------------------------------------------------------------
% Schaefer 100
% omega=0.1
% 0-back         % gamma=1      % gamma=1.1       % gamma=1.2        % gamma=1.3     % gamma=1.4       % gamma=1.5       % gamma=1.6     %         
%c=[2,5,10,11];  % r=*          % r=#             % r=#              % r=#           % r=#             % r=#             % r=#           %  
%c=[2,10,5,11];  % r=*          % r=#             % r=#              % r=#           % r=#             % r=#             % r=#           %    
%c=[2,11,5,10];  % r=*          % r=#             % r=#              % r=#           % r=#             % r=#             % r=#           % 
                 % r_mean=*     % r_mean=#        % r_mean=#         % r_mean=#      % r_mean=#        % r_mean=#        % r_mean=#      %   

%--------------------------------------------------------------------------
% Schaefer 100
% omega=0.5
% 2-back         % gamma=1    % %       
%c=[1,4,7,8];    % r=#        % %    
%c=[1,7,4,8];    % r=#         % % 
%c=[1,8,4,7];    % r=#         % %        
                 % r_mean=#    % % 
%--------------------------------------------------------------------------
% Schaefer 100
% omega=0.5
% 0-back         % gamma=1      %       
%c=[2,5,10,11];  % r=#          % 
%c=[2,10,5,11];  % r=#          %    
%c=[2,11,5,10];  % r=#          %  
                 % r_mean=#     %  


%--------------------------------------------------------------------------
% Kong 200
% omega=0.01       
% 2-back         % gamma=1    % gamma=1.5        % gamma=2      
%c=[1,4,7,8];    % r=#        % r=0.3072         % r=0.2850   
%c=[1,7,4,8];    % r=#        % r=0.2593         % r=0.2875
%c=[1,8,4,7];    % r=#        % r=0.3477         % r=0.4070
                              % r_mean=0.3047    % r_mean=0.3265
%--------------------------------------------------------------------------
% Kong 200
% omega=0.01
% 0-back         % gamma=1      % gamma=1.5      % gamma=2
%c=[2,5,10,11];  % r=#          % r=0.2624        % r=0.3159
%c=[2,10,5,11];  % r=#          % r=0.4256         % r=0.3718
c=[2,11,5,10];  % r=#          % r=0.3266          % r=0.3256
                                % r_mean=0.3382    % r_mean=0.3378


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
half_1
half_2
r=corr(half_1,half_2)
% 
% data_path = fileparts(mfilename('fullpath'));
% NMI_indi_path=fullfile(data_path,['Results/real_modularity/',num2str(gamma),'/NMI_indi_real_modularity']);
% save(NMI_indi_path,'NMI_indi_real_modularity')


