% Test CommuDetectLBM.m
clear
clc
close all
load('/Users/lingbinbian/Documents/MCD_revision_v2/MCD1.0/CommuDetect/Results/synthetic_LBM/DIIV10/n0.3162/grouplevel_data.mat')


A=zeros(200,200);
B=zeros(200,200);
C=zeros(200,200);
D=zeros(200,200);

A(1:100,1:100)=group_adj{1,1}{1,1};
A(1:100,101:200)=group_adj{2,1}{1,1};
A(101:200,1:100)=group_adj{3,1}{1,1};
A(101:200,101:200)=group_adj{4,1}{1,1};

A(1:100,1:100)=group_adj{5,1}{1,1};
A(1:100,101:200)=group_adj{6,1}{1,1};
A(101:200,1:100)=group_adj{7,1}{1,1};
A(101:200,101:200)=group_adj{8,1}{1,1};

C(1:100,1:100)=group_adj{9,1}{1,1};
C(1:100,101:200)=group_adj{10,1}{1,1};
C(101:200,1:100)=group_adj{11,1}{1,1};
C(101:200,101:200)=group_adj{12,1}{1,1};

C(1:100,1:100)=group_adj{13,1}{1,1};
C(1:100,101:200)=group_adj{14,1}{1,1};
C(101:200,1:100)=group_adj{15,1}{1,1};
C(101:200,101:200)=group_adj{16,1}{1,1};


label_vector=CommuDetectLBM([A,B;C,D]);

