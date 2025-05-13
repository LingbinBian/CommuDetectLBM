function [colormap_type] = color_type(v)
% Colormap
% v: a vector of labels

color_1=[0.69,0.19,0.38];
color_2=[1.00,0.89,0.52];
color_3=[0.25,0.41,0.88];
color_4=[0.00,0.79,0.34];
color_5=[0.63,0.32,0.18];
color_6=[0.50,0.54,0.53];
color_7=[0.53,0.81,0.92];
color_8=[0.63,0.40,0.83];
color_9=[0.74,0.56,0.56];
color_10=[0.16,0.14,0.13];
color_11=[1.00,0.50,0.31];
color_12=[1.00,0.00,0.00];
color_13=[0.39,0.19,0.38];
color_14=[0.50,0.89,0.52];
color_15=[0.12,0.41,0.88];
color_16=[0.31,0.32,0.18];
color_17=[1.00,0.00,1.00];
color_18=[0.25,0.54,0.53];
color_19=[0.27,0.21,0.92];
color_20=[0.17,0.57,0.17];
color_map_RGB=[color_1;color_2;color_3;color_4;color_5;color_6;color_7;color_8;color_9;color_10;color_11;color_12;color_13;color_14;color_15;color_16;color_17;color_18;color_19;color_20];
colormap_type=color_map_RGB(unique(v),:);

end






