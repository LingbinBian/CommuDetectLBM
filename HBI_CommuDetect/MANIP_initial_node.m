% Initialize .node 

clear
clc

centroid_Schaefer100=csvread('Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm.Centroid_RAS.csv',1,2);

spar_level=10;
Node_ROI=dlmread(['Hierarchical/real/LR/sparsity_',num2str(spar_level),'/Node_ROI_t',num2str(41),'.node']);
Node_ROI(:,1:3)=centroid_Schaefer100;
Node_ROI(:,6)=1:100;

dlmwrite(['Hierarchical/real/LR/sparsity_',num2str(spar_level),'/Node_ROI_t',num2str(41),'.node'],Node_ROI,'delimiter','\t')
