% Generate some sample data
x = randn(100,1);
y = x + randn(100,1)*0.5; % y is correlated with x plus some noise

% Create scatter plot
% scatter(x, y);
% xlabel('X Variable');
% ylabel('Y Variable');
% title('Scatter Plot Showing Correlation');
% grid on;

% Using the first example's data
scatter(x, y,12,"filled");
hold on;

% Calculate and plot regression line
p = polyfit(x, y, 1);
yfit = polyval(p, x);
plot(x, yfit, 'k-', 'LineWidth', 1);

% Add correlation coefficient to plot
r = corrcoef(x, y);
r = r(1,2);
text(min(x)+1, 0.9*max(y), sprintf('r = %.2f', r), ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left','fontsize', 16, 'fontname', 'times');



xlabel('X Variable','fontsize', 16, 'fontname', 'times');
ylabel('Y Variable','fontsize', 16, 'fontname', 'times');
title('Correlation','fontsize', 16, 'fontname', 'times');
set(gca, 'linewidth', 1.2, 'fontsize', 12, 'fontname', 'times')
set(gca,'box','on')
set(gcf,'unit','normalized','position',[0.3,0.22,0.115,0.185]);
grid on;
hold off;