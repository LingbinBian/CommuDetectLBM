function[]=visual_labels(labels)
% This function visualizes the latent labels.
    imagesc(labels)    
    xlabel('Label','fontsize',16)
    ylabel('Node','fontsize',16)
    % colormap(parula(max(K_min)));
   %  colormap(lines(max(K_min)));
    colormap(color_type(labels));
    colorbar_community(max(max(labels)))
   % colorbar_community_K(labels)
    set(gca,'xtick',[])
    set(gca, 'linewidth', 1.2, 'fontsize', 16, 'fontname', 'times')
end