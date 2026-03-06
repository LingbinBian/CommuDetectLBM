function [] = correlation_plot(x,x_name,y,y_name)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% Generate some sample data
% x = randn(100,1);
% y = x + randn(100,1)*0.5; % y is correlated with x plus some noise

% Create scatter plot
% scatter(x, y);
% xlabel('X Variable');
% ylabel('Y Variable');
% title('Scatter Plot Showing Correlation');
% grid on;

% Using the first example's data
scatter(x, y, 8, "filled");
hold on;

% Calculate and plot regression line
p = polyfit(x, y, 1);
yfit = polyval(p, x);
plot(x, yfit, 'k-', 'LineWidth', 1);

% Add correlation coefficient to plot
r = corrcoef(x, y);
r = r(1,2);
text(min(x)+0.1*min(x), 0.9*max(y), sprintf('r = %.2f', r), ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left','fontsize', 16, 'fontname', 'times');



xlabel(x_name,'fontsize', 16, 'fontname', 'times');
ylabel(y_name,'fontsize', 16, 'fontname', 'times');
title('Correlation','fontsize', 16, 'fontname', 'times');
set(gca, 'linewidth', 1.2, 'fontsize', 16, 'fontname', 'times')
set(gca,'box','on')
set(gcf,'unit','normalized','position',[0.3,0.22,0.12,0.18]);
%grid on;
hold off;
end