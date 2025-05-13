% This script saves all the figures of community labels to the destination
% file Pictures_community_labels
% Version 1.0
% 1-May-2025
% Copyright (c) 2025, Lingbin Bian

clear
clc
close all

% destination file

pic_path='../Pictures_community_labels/synthetic/LBM/DIIV30/';
  
% save hierarchical synthetic
cd 'Results/synthetic_LBM/DIIV30/n0.3162';
open('assign_prob_max_state_1.fig')
saveas(gcf,['../../../../',pic_path,'n0.3162/assign_prob_max_state_1.svg'])
close all

open('assign_prob_max_state_2.fig')
saveas(gcf,['../../../../',pic_path,'n0.3162/assign_prob_max_state_2.svg'])
close all

open('assign_prob_max_state_3.fig')
saveas(gcf,['../../../../',pic_path,'n0.3162/assign_prob_max_state_3.svg'])
close all


open('assign_prob_state_1.fig')
saveas(gcf,['../../../../',pic_path,'n0.3162/assign_prob_state_1.svg'])
close all

open('assign_prob_state_2.fig')
saveas(gcf,['../../../../',pic_path,'n0.3162/assign_prob_state_2.svg'])
close all

open('assign_prob_state_3.fig')
saveas(gcf,['../../../../',pic_path,'n0.3162/assign_prob_state_3.svg'])
close all



open('group_latent_labels_1.fig')
saveas(gcf,['../../../../',pic_path,'n0.3162/group_latent_labels_1.svg'])
close all

open('group_latent_labels_2.fig')
saveas(gcf,['../../../../',pic_path,'n0.3162/group_latent_labels_2.svg'])
close all

open('group_latent_labels_3.fig')
saveas(gcf,['../../../../',pic_path,'n0.3162/group_latent_labels_3.svg'])
close all


open('Labels_esti_1.fig')
saveas(gcf,['../../../../',pic_path,'n0.3162/Labels_esti_1.svg'])
close all

open('Labels_esti_2.fig')
saveas(gcf,['../../../../',pic_path,'n0.3162/Labels_esti_2.svg'])
close all

open('Labels_esti_3.fig')
saveas(gcf,['../../../../',pic_path,'n0.3162/Labels_esti_3.svg'])
close all

open('Labels_switched_1.fig')
saveas(gcf,['../../../../',pic_path,'n0.3162/Labels_switched_1.svg'])
close all

open('Labels_switched_2.fig')
saveas(gcf,['../../../../',pic_path,'n0.3162/Labels_switched_2.svg'])
close all

open('Labels_switched_3.fig')
saveas(gcf,['../../../../',pic_path,'n0.3162/Labels_switched_3.svg'])
close all

open('Labels_true_1.fig')
saveas(gcf,['../../../../',pic_path,'n0.3162/Labels_true_1.svg'])
close all

open('Labels_true_2.fig')
saveas(gcf,['../../../../',pic_path,'n0.3162/Labels_true_2.svg'])
close all

open('Labels_true_3.fig')
saveas(gcf,['../../../../',pic_path,'n0.3162/Labels_true_3.svg'])
close all


