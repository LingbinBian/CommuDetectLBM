function []=colorbar_community_K(z_uni)
% Define the colorbar for communities.
%
% Version 1.0 
% 06-Nov-2019
% Copyright (c) 2019, Lingbin Bian

colorbar('Ticks',z_uni,...
         'TickLabels',{num2str(z_uni)},...
         'fontsize',16);

end