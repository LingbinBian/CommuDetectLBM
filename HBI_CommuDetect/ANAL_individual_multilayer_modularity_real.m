% Community detection using multilayer modularity (real fMRI data)
% individual-level analysis
%
% Version 1.0
% 23-Feb-2026
% Copyright (c) 2026, Lingbin Bian
% -------------------------------------------------------------------------
clear
clc
close all
% -------------------------------------------------------------------------
% Load data
datatype=1;
n_c=1;

switch n_c
    case 1               % 2-back 
        c=[1,4,7,8];     % Tool & Face; Body & Place 
    case 2
        c=[1,7,4,8];     % Tool & Body; Face & Place 
    case 3
        c=[1,8,4,7];     % Tool & Place; Face & Body 
    case 4               % 0-back 
        c=[2,5,10,11];   % Tool & Face; Body & Place
    case 5
        c=[2,10,5,11];   % Tool & Body; Face & Place 
    case 6
        c=[2,11,5,10];   % Tool & Place; Face & Body 
end



subjects=load('subject.txt');
atlas=2; % 1: Scheafer's atlas, 2: Kong's atlas
session_n=1;
n_s='';
vari='';
hrf_ind='';
if session_n==1
    state_t=[41,76,107,140,175,206,239,278,306,334,375];
elseif session_n==2
    state_t=[49,77,99,139,175,209,236,275,305,334,376];
end
K_state=[6,6,6,6,6,6,6,6,6,6,6];  % make an assumption of K
W=10;

load('Results/real_Bayes_opt_MM/opt_para_multi_modularity.mat')
para_M=opt_para_MM;  % gamma

gamma=para_M{1,n_c};
omega=para_M{2,n_c};


N_subj=100;

N_state=length(state_t);

group_adj=cell(N_subj,1); % group of adjacency matrix
Q=zeros(N_subj,N_state);
z_g=cell(N_subj,N_state); % group of estimated z

% K_chain=cell(N_subj,N_state);
% Z_chain=cell(N_subj,N_state);

for s=1:N_subj
    fprintf('Adjacency matrix of subject: %d\n',s)
    subid=num2str(subjects(s));
    [group_adj{s,1},true_latent,true_latent_sub]=local_adj(datatype,atlas,subid,session_n,n_s,state_t,K_state,W,vari,hrf_ind);
end

% -------------------------------------------------------------------------
% Adjacency matrix


M=10; % margin size
if atlas==1
    N=100; 
else
    N=200; 
end
group_adj_3D=cell(N_state,N_subj);
communities=zeros(N,N_subj,N_state);
Q_multi=zeros(1,N_state);

for s=1:N_subj
    fprintf('subject: %d\n',s)
    for ds=1:N_state
       group_adj_3D{ds,s}=group_adj{s,1}{1,ds};
    end
end

for ds=1:N_state
    [communities(:,:,ds)] = multilayer_modularity(group_adj_3D(ds,:),3.92,0.0379);
end

for s=1:N_subj
    fprintf('subject: %d\n',s)
    for ds=1:N_state
        %[Q_multi(1,ds), communities(:,:,ds)] = multilayer_modularity(group_adj_3D, r, w, 'categorical');
        %[z_g{s,ds},Q(s,ds)] = modularity_und(group_adj{s,1}{1,ds},r);   
        z_g{s,ds}=communities(:,s,ds);       
        
    end
end

% normalized mutual information

half_1=zeros(N_subj,1);
half_2=zeros(N_subj,1);

for s=1:N_subj
        half_1(s,1)=nmi(z_g{s,c(1)},z_g{s,c(2)});
        half_2(s,1)=nmi(z_g{s,c(3)},z_g{s,c(4)});
end

r=corr(half_1,half_2)

% Schaefer's atlas
% Mulilayer modularity
% c1=0.4309     gamma=2.5940   omega=0.0392
% c2=
% c3=
% c4=
% c5=
% c6=











