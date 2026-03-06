% Community detection using LBM (real fMRI data)
% individual-level analysis
%
% Version 1.0
% 25-Feb-2026
% Copyright (c) 2026, Lingbin Bian
% -------------------------------------------------------------------------
clear
clc
close all
% -------------------------------------------------------------------------
% Load data
datatype=1;

atlas=1; % 1: Scheafer's atlas, 2: Kong's atlas

c_comb=1;

switch c_comb
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

N_sub=100;
subjects=load('subject.txt');

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
% 
% if atlas==1
%     load(['Results/real_Bayes_opt_LBM/combination',num2str(c_comb),'.mat'])
%     K_max,nu,rho,kappa_sq
% else
%     load(['Results/real_Kong_Bayes_opt_LBM/combination',num2str(c_comb),'.mat'])
%     gamma_M=bestGamma;  % gamma
% end


K_max=20;
nu=3;
rho=0.02;
kappa_sq=1;

N_subj=100;

N_state=length(state_t);

group_adj=cell(N_subj,1); % group of adjacency matrix
Q=zeros(N_subj,N_state);
z_g=cell(N_subj,4); % group of estimated z

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

for s=1:N_subj
    fprintf('subject: %d\n',s)
    for ds=1:4
      [z_g{s,ds},Q(s,ds)] = CommuDetectLBM(group_adj{s,1}{1,ds},K_max,nu,rho,kappa_sq);   
      %[z_g{s,ds},Q(s,ds)] = CommuDetectLBM(group_adj{s,1}{1,c(ds)},20,3,0.02,1);  
        
      
    end
end

% normalized mutual information

half_1=zeros(N_sub,1);
half_2=zeros(N_sub,1);

for s=1:N_sub
        half_1(s,1)=nmi(z_g{s,1},z_g{s,2});
        half_2(s,1)=nmi(z_g{s,3},z_g{s,4});
end

r=corr(half_1,half_2)

if atlas==1
    data_path = fileparts(mfilename('fullpath'));
    para_path=fullfile(data_path,['Results/real_Bayes_opt_LBM','/split_corre']);
    save(para_path,'nu','rho','kappa_sq','K_max','r')
else
    data_path = fileparts(mfilename('fullpath')); 
    para_path=fullfile(data_path,['Results/real_Kong_Bayes_opt_LBM','/split_corre']);
    save(para_path,'nu','rho','kappa_sq','K_max','r')
end

% Schaefer's atlas
% Modularity
% c1=0.3845    gamma=2.8450
% c2=0.3870    gamma=2.8472
% c3=0.4508    gamma=2.8927
% c4=0.3726    gamma=3.2283
% c5=0.4560    gamma=1.5021
% c6=0.3753    gamma=3.4013

% Kong's atlas
% Modularity
% c1=            
% c2=0.3937          gamma=2.8714









