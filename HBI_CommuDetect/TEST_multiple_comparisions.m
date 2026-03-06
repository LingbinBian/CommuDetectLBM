% Set the random seed for reproducibility
rng(0);

% Number of subjects
numSubjects = 100;

% Simulate NMI values
nmi_latent_block = rand(numSubjects, 1); % Latent Block Model values
nmi_modularity_gamma1 = rand(numSubjects, 1) * 0.75; 
nmi_modularity_gamma1_5 = rand(numSubjects, 1) * 0.6;
nmi_modularity_gamma2 = rand(numSubjects, 1) * 0.5;

% Simulating multilayer modularity values
nmi_multilayer_1_0_01 = rand(numSubjects, 1) * 0.7;
nmi_multilayer_1_0_1 = rand(numSubjects, 1) * 0.65;
nmi_multilayer_1_0_5 = rand(numSubjects, 1) * 0.55;

% Combine NMI data for modularity into a single array
modularity_nmi = [
    nmi_modularity_gamma1; 
    nmi_modularity_gamma1_5; 
    nmi_modularity_gamma2;
    nmi_multilayer_1_0_01; 
    nmi_multilayer_1_0_1; 
    nmi_multilayer_1_0_5
];

% Combine labels
modularity_labels = [
    repmat({'Modularity (γ=1)'}, numSubjects, 1); 
    repmat({'Modularity (γ=1.5)'}, numSubjects, 1); 
    repmat({'Modularity (γ=2)'}, numSubjects, 1);
    repmat({'Multilayer Modularity (ω=0.01)'}, numSubjects, 1);
    repmat({'Multilayer Modularity (ω=0.1)'}, numSubjects, 1);
    repmat({'Multilayer Modularity (ω=0.5)'}, numSubjects, 1)
];

% Combine all for the boxplot
all_data = [nmi_latent_block; modularity_nmi];
all_labels = [repmat({'Latent Block'}, numSubjects, 1); modularity_labels];

% Visualization
figure;

% Boxplot
boxplot(all_data, all_labels, 'Colors', 'k', 'BoxStyle', 'outline', 'Whisker', 1.5);
title('NMI Comparison: Latent Block vs Modularity');
xlabel('Method');
ylabel('NMI Value');

% Set significance level indicators
hold on;
text(1.2, 0.7, '***', 'FontSize', 12, 'Color', 'r', 'HorizontalAlignment', 'center');
text(1.8, 0.65, '***', 'FontSize', 12, 'Color', 'r', 'HorizontalAlignment', 'center');
text(1.2, 0.5, '*', 'FontSize', 12, 'Color', 'r', 'HorizontalAlignment', 'center');

% Customize axes
ylim([0 1]);
grid on;
set(gca, 'FontSize', 12, 'XTickLabelRotation', 45);

% Create a Legend
legend({'Latent Block Model', ...
        'Modularity (γ=1)', ...
        'Modularity (γ=1.5)', ...
        'Modularity (γ=2)', ...
        'Multilayer Modularity (ω=0.01)', ...
        'Multilayer Modularity (ω=0.1)', ...
        'Multilayer Modularity (ω=0.5)'}, ...
        'Location', 'northeastoutside', 'Orientation', 'horizontal');

hold off;