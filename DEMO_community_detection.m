% Demonstration of CommuDetectLBM.m and CommuDetectGroup
% Version 1.0
% 10-April-2025
% Copyright (c) 2025, Lingbin Bian, 
% If the users have any questions, please contact via lingbin.bian@gmail.com
% -------------------------------------------------------------------------
clear
clc

% examples of adjacency matrices
A1=[9,9,0,0,0,0,0,0,0;
   9,9,0,0,0,0,0,0,0;
   0,0,8,8,8,0,0,0,0;
   0,0,8,8,8,0,0,0,0;
   0,0,8,8,8,0,0,0,0;
   0,0,0,0,0,7,7,7,7;
   0,0,0,0,0,7,7,7,7;
   0,0,0,0,0,7,7,7,7;
   0,0,0,0,0,7,7,7,7;];

A2=[9,9,9,9,9,1,1,1,1;
   9,9,9,9,9,1,1,1,1;
   9,9,9,9,9,1,1,1,1;
   9,9,9,9,9,1,1,1,1;
   9,9,9,9,9,1,1,1,1;
   1,1,1,1,1,7,7,7,7;
   1,1,1,1,1,7,7,7,7;
   1,1,1,1,1,7,7,7,7;
   1,1,1,1,1,7,7,7,7;];

A3=[7,7,0,0,0,0,0,0,0;
    7,7,0,0,0,0,0,0,0;
    0,0,25,25,0,0,0,0,0;
    0,0,25,25,0,0,0,0,0;
    0,0,0,0,19,19,19,0,0;
    0,0,0,0,19,19,19,0,0;
    0,0,0,0,19,19,19,0,0;
    0,0,0,0,0,0,0,25,25;
    0,0,0,0,0,0,0,25,25;];

% Individual-level community detection
z1=CommuDetectLBM(A1);
z2=CommuDetectLBM(A2);
z3=CommuDetectLBM(A3);

% Group-level community detection
Z=[z1,z2,z3];
z_group=CommuDetectGroup(Z);


