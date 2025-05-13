% Hierarchical Bayesian inference (HBI)
% group-level anlaysis: modelling connectivity
%
% Version 1.0
% 30-April-2025
% Copyright (c) 2025, Lingbin Bian

clear
clc
close all



datatype=0;

if datatype==0
  % signal to noise ratio SNR=10log10(sigma_signal^2/sigma_noise^2);
  % sigma_signal=1, sigma_noise=n_s
   N=100;
   n_s=0.3162;  % 10dB
  % n_s=0.5623;  % 5dB
  % n_s=1;       % 0dB
  % n_s=1.7783;  % -5dB
  % n_s=3.1623;  % -10dB
  % n_s=5.6234;  % -15dB
  % n_s=10;     % -20dB
  % n_s=17.7828; % -25dB

    DIIV=20;
    load(['Results/synthetic_LBM/DIIV',num2str(DIIV),'/n',num2str(n_s),'/grouplevel_data.mat']);

elseif datatype==1
    atlas=2;
    if atlas==1
        N=100;
    else
        N=200;
    end
    session_n=1;
    if atlas==1

        if session_n==1
            load('Results/real_LBM/LR/grouplevel_data.mat');    
        elseif session_n==2
            load('Results/real_LBM/RL/grouplevel_data.mat'); 
        end
    else
        if session_n==1
            load('Results/real_Kong_LBM/LR/grouplevel_data.mat');    
        elseif session_n==2
            load('Results/real_Kong_BM/RL/grouplevel_data.mat'); 
        end
    end
end

% Modelling connectivity
[esti_connectmean,esti_connectvariance]=connectivity_modelling(group_adj,N,state_t,20);

% point estimate
ave_adj=cell(1,N_state);  % averaged adjacency matrix
variance_adj=cell(1,N_state);  % variance matrix

adj_mean=zeros(N_subj,1);
for t=1:N_state
    for i=1:N
        for j=1:N
            for s=1:N_subj
                adj_mean(s,1)=group_adj{s,1}{1,t}(i,j);
            end
            ave_adj{1,t}(i,j)=mean(adj_mean);
        end
    end
end

adj_vari=zeros(N_subj,1);
for t=1:N_state
    for i=1:N
        for j=1:N
            for s=1:N_subj
                adj_vari(s,1)=group_adj{s,1}{1,t}(i,j);
            end
            variance_adj{1,t}(i,j)=var(adj_vari);
        end
    end
end


% Visualize connectivity mean ---------------------------------------------
figure

for s=1:N_state
   if datatype==0  
      subplot(1,3,s)
   elseif datatype==1
      subplot(3,4,s)
   end 
    imagesc(esti_connectmean{s})
    colormap(turbo);
    colorbar
    title(['Post mean, State',' ',num2str(s)],'fontsize',16) 
    xlabel('Node','fontsize',16) 
    ylabel('Node','fontsize',16)
      
    set(gca, 'linewidth', 1.2, 'fontsize', 16, 'fontname', 'times')
end

if datatype==0
    set(gcf,'unit','normalized','position',[0.3,0.2,0.5,0.18]);
elseif datatype==1
    set(gcf,'unit','normalized','position',[0.3,0.2,0.7,0.6]);
end
data_path = fileparts(mfilename('fullpath'));       
if datatype==0
    saveas(gcf,['Results/synthetic_LBM/DIIV',num2str(DIIV),'/n',num2str(n_s),'/group_connectmean_estimation.fig'])
    saveas(gcf,['Results/synthetic_LBM/DIIV',num2str(DIIV),'/n',num2str(n_s),'/group_connectmean_estimation.svg'])
elseif datatype==1
    if atlas==1
        if session_n==1
            L_results_path=fullfile(data_path,'Results/real_LBM/LR/LR_meanconnectivity');
            save(L_results_path,'esti_connectmean','ave_adj');
            saveas(gcf,['Results/real_LBM/LR/','group_connectmean_estimation.fig'])
        elseif session_n==2
            R_results_path=fullfile(data_path,'Results/real_LBM/RL/RL_meanconnectivity');
            save(R_results_path,'esti_connectmean','ave_adj');
            saveas(gcf,['Results/real_LBM/RL/','group_connectmean_estimation.fig'])
        end
    else
        if session_n==1
            L_results_path=fullfile(data_path,'Results/real_Kong_LBM/LR/LR_meanconnectivity');
            save(L_results_path,'esti_connectmean','ave_adj');
            saveas(gcf,['Results/real_Kong_LBM/LR/','group_connectmean_estimation.fig'])
        elseif session_n==2
            R_results_path=fullfile(data_path,'Results/real_Kong_LBM/RL/RL_meanconnectivity');
            save(R_results_path,'esti_connectmean','ave_adj');
            saveas(gcf,['Results/real_Kong_LBM/RL/','group_connectmean_estimation.fig'])
        end
    end
end

% Visualize connectivity variance -----------------------------------------
figure

for s=1:N_state
   if datatype==0  
      subplot(1,3,s)
   elseif datatype==1
      subplot(3,4,s)
   end 
    imagesc(esti_connectvariance{s})
    colormap(parula);
    colorbar
    title(['Post variance, State',' ',num2str(s)],'fontsize',16) 
    xlabel('Node','fontsize',16) 
    ylabel('Node','fontsize',16)
      
    set(gca, 'linewidth', 1.2, 'fontsize', 16, 'fontname', 'times')
end

if datatype==0
    set(gcf,'unit','normalized','position',[0.3,0.2,0.5,0.18]);
elseif datatype==1
    set(gcf,'unit','normalized','position',[0.3,0.2,0.7,0.6]);
end
        
if datatype==0
    saveas(gcf,['Results/synthetic_LBM/DIIV',num2str(DIIV),'/n',num2str(n_s),'/group_connectvariance_estimation.fig'])
    saveas(gcf,['Results/synthetic_LBM/DIIV',num2str(DIIV),'/n',num2str(n_s),'/group_connectvariance_estimation.svg'])
elseif datatype==1
    if atlas==1
        if session_n==1
            saveas(gcf,['Results/real_LBM/LR/','group_connectvariance_estimation.fig'])
        elseif session_n==2
            saveas(gcf,['Results/real_LBM/RL/','group_connectvariance_estimation.fig'])
        end
    else
        if session_n==1
            saveas(gcf,['Results/real_Kong_LBM/LR/','group_connectvariance_estimation.fig'])
        elseif session_n==2
            saveas(gcf,['Results/real_Kong_LBM/RL/','group_connectvariance_estimation.fig'])
        end
        
    end
end

% Visualize mean matrix (point estimate)-----------------------------------
figure

for s=1:N_state
   if datatype==0  
      subplot(1,3,s)
   elseif datatype==1
      subplot(3,4,s)
   end 
    imagesc(ave_adj{s})
    colormap(turbo);
    colorbar
    title(['Empirical mean, State',' ',num2str(s)],'fontsize',16) 
    xlabel('Node','fontsize',16) 
    ylabel('Node','fontsize',16)
      
    set(gca, 'linewidth', 1.2, 'fontsize', 16, 'fontname', 'times')
end

if datatype==0
    set(gcf,'unit','normalized','position',[0.3,0.2,0.5,0.18]);
elseif datatype==1
    set(gcf,'unit','normalized','position',[0.3,0.2,0.7,0.6]);
end
        
if datatype==0
    saveas(gcf,['Results/synthetic_LBM/DIIV',num2str(DIIV),'/n',num2str(n_s),'/group_ave_adj.fig'])
    saveas(gcf,['Results/synthetic_LBM/DIIV',num2str(DIIV),'/n',num2str(n_s),'/group_ave_adj.svg'])
elseif datatype==1
    if atlas==1
        if session_n==1
            saveas(gcf,['Results/real_LBM/LR/','group_ave_adj.fig'])
        elseif session_n==2
            saveas(gcf,['Results/real_LBM/RL/','group_ave_adj.fig'])
        end
    else
        if session_n==1
            saveas(gcf,['Results/real_Kong_LBM/LR/','group_ave_adj.fig'])
        elseif session_n==2
            saveas(gcf,['Results/real_Kong_LBM/RL/','group_ave_adj.fig'])
        end
    end
end

% Visualize variance matrix (point estimate) ------------------------------
figure

for s=1:N_state
   if datatype==0  
      subplot(1,3,s)
   elseif datatype==1
      subplot(3,4,s)
   end 
    imagesc(variance_adj{s})
    colormap(parula);
    colorbar
    title(['Empirical variance, State',' ',num2str(s)],'fontsize',16) 
    xlabel('Node','fontsize',16) 
    ylabel('Node','fontsize',16)
      
    set(gca, 'linewidth', 1.2, 'fontsize', 16, 'fontname', 'times')
end

if datatype==0
    set(gcf,'unit','normalized','position',[0.3,0.2,0.5,0.18]);
elseif datatype==1
    set(gcf,'unit','normalized','position',[0.3,0.2,0.7,0.6]);
end
        
if datatype==0
    saveas(gcf,['Results/synthetic_LBM/DIIV',num2str(DIIV),'/n',num2str(n_s),'/group_variance_adj.fig'])
    saveas(gcf,['Results/synthetic_LBM/DIIV',num2str(DIIV),'/n',num2str(n_s),'/group_variance_adj.svg'])
elseif datatype==1
    if atlas==1
        if session_n==1
            saveas(gcf,['Results/real_LBM/LR/','group_variance_adj.fig'])
        elseif session_n==2
            saveas(gcf,['Results/real_LBM/RL/','group_variance_adj.fig'])
        end
    else
        if session_n==1
            saveas(gcf,['Results/real_Kong_LBM/LR/','group_variance_adj.fig'])
        elseif session_n==2
            saveas(gcf,['Results/real_Kong_LBM/RL/','group_variance_adj.fig'])
        end
    end
end
    
% -------------------------------------------------------------------------
% % comparison
connect_mean_row=zeros(N_state,N*N-N);
ave_adj_row=zeros(N_state,N*N-N);


for s=1:N_state
    count=0;
    figure
    for i=1:N
        for j=1:N
            if i~=j
               % scatter(esti_connectmean{s}(i,j),ave_adj{s}(i,j),8,'filled')
               % hold on
               connect_mean_row(s,count+1)=esti_connectmean{s}(i,j);
               ave_adj_row(s,count+1)=ave_adj{s}(i,j);
               count=count+1;

            end
        end
    end
    count=0;
    % R_m=corrcoef(connect_mean_row(s,:),ave_adj_row(s,:));
    % R_mean(s,1)=R_m(1,2);
    % 
    % if datatype==0
    %    x=linspace(-0.1,1);
    %    y=linspace(-0.1,1);
    %    plot(x,y,'k--','linewidth', 1.2);
    % elseif datatype==1
    %    x=linspace(-0.1,0.6);
    %    y=linspace(-0.1,0.6);
    %    plot(x,y,'k--','linewidth', 1.2);
    % end
    % tableData_mean = table(connect_mean_row(s,:)', ave_adj_row(s,:)', ...
    % 'VariableNames', {'Post', 'Real'});
    % 
    % [R,PValue] = corrplot(tableData_mean);
    % hold on
    % set(gca, 'linewidth', 1.2, 'fontsize', 12, 'fontname', 'times')
    % set(gca,'box','on')
    % 
    % title(['Mean, State ',num2str(s)],'fontsize',12) 
 
    correlation_plot(connect_mean_row(s,:),'Posterior mean',ave_adj_row(s,:),'Empirical mean');



    
    if datatype==0
        saveas(gcf,['Results/synthetic_LBM/DIIV',num2str(DIIV),'/n',num2str(n_s),'/comparison_mean_',num2str(s),'.fig'])
        saveas(gcf,['Results/synthetic_LBM/DIIV',num2str(DIIV),'/n',num2str(n_s),'/comparison_mean_',num2str(s),'.svg'])
    elseif datatype==1
        if atlas==1
            if session_n==1
                saveas(gcf,['Results/real_LBM/LR/','comparison_mean.fig'])
            elseif session_n==2
                saveas(gcf,['Results/real_LBM/RL/','comparison_mean.fig'])
            end
        else
            if session_n==1
                saveas(gcf,['Results/real_Kong_LBM/LR/','comparison_mean.fig'])
            elseif session_n==2
                saveas(gcf,['Results/real_Kong_LBM/RL/','comparison_mean.fig'])
            end
        end
    end
end


% variance

connect_variance_row=zeros(N_state,N*N-N);
variance_adj_row=zeros(N_state,N*N-N);

for s=1:N_state
    count=0;
    figure
    for i=1:N
        for j=1:N
            if i~=j          
               connect_variance_row(s,count+1)=esti_connectvariance{s}(i,j);
               variance_adj_row(s,count+1)=variance_adj{s}(i,j);
               count=count+1;

            end
        end
    end
    count=0;

    % tableData_variance = table(connect_variance_row(s,:)', variance_adj_row(s,:)', ...
    % 'VariableNames', {'Post', 'Real'});
    % 
    % [R,PValue] = corrplot(tableData_variance);
    % hold on
    % set(gca, 'linewidth', 1.2, 'fontsize', 12, 'fontname', 'times')
    % set(gca,'box','on')
    % 
    % title(['Variance, State ',num2str(s)],'fontsize',12) 
    correlation_plot(connect_variance_row(s,:),'Posterior variance',variance_adj_row(s,:),'Empirical variance');
 



    
    if datatype==0
        saveas(gcf,['Results/synthetic_LBM/DIIV',num2str(DIIV),'/n',num2str(n_s),'/comparison_variance_',num2str(s),'.fig'])
        saveas(gcf,['Results/synthetic_LBM/DIIV',num2str(DIIV),'/n',num2str(n_s),'/comparison_variance_',num2str(s),'.svg'])
    elseif datatype==1
        if atlas==1
            if session_n==1
                saveas(gcf,['Results/real_LBM/LR/','comparison_variance.fig'])
            elseif session_n==2
                saveas(gcf,['Results/real_LBM/RL/','comparison_variance.fig'])
            end
        else
            if session_n==1
                saveas(gcf,['Results/real_Kong_LBM/LR/','comparison_variance.fig'])
            elseif session_n==2
                saveas(gcf,['Results/real_Kong_LBM/RL/','comparison_variance.fig'])
            end
            
        end
    end
end

% connect_variance_row=zeros(N_state,N*N-N);
% variance_adj_row=zeros(N_state,N*N-N);
% R_variance=zeros(N_state,1);
% 
% for s=1:N_state
%     count=0;
%     figure
%     for i=1:N
%         for j=1:N
%             if i~=j
%                scatter(esti_connectvariance{s}(i,j),variance_adj{s}(i,j),8,'filled')
%                hold on
%                connect_variance_row(s,count+1)=esti_connectvariance{s}(i,j);
%                variance_adj_row(s,count+1)=variance_adj{s}(i,j);
%                count=count+1;
%             end
%         end
%     end
%     count=0;
%     R_v=corrcoef(connect_variance_row(s,:),variance_adj_row(s,:));
%     R_variance(s,1)=R_v(1,2);
% 
%     if datatype==0
%        x=linspace(-0.01,0.1);
%        y=linspace(-0.01,0.1);
%        plot(x,y,'k--','linewidth', 1.2);
%     elseif datatype==1
%        x=linspace(-0.01,0.22);
%        y=linspace(-0.01,0.22);
%        plot(x,y,'k--','linewidth', 1.2);
%     end
%     hold on
% 
% 
%     set(gca, 'linewidth', 1.2, 'fontsize', 16, 'fontname', 'times')
%     set(gca,'box','on')
%     title(['State ',num2str(s)],'fontsize',16) 
%     xlabel('Posterior sample','fontsize',16) 
%     ylabel('Data variance','fontsize',16)
% 
% 
% if datatype==0
%     set(gcf,'unit','normalized','position',[0.3,0.22,0.115,0.185]);
% elseif datatype==1
%     set(gcf,'unit','normalized','position',[0.3,0.22,0.115,0.185]);
% end
% 
% if datatype==0
%     saveas(gcf,['Results/synthetic/','n',num2str(n_s),'/comparison_variance_',num2str(s),'.fig'])
% elseif datatype==1
%     if session_n==1
%         saveas(gcf,['Results/real/LR/','comparison_variance.fig'])
%     elseif session_n==2
%         saveas(gcf,['Results/real/RL/','comparison_variance.fig'])
%     end
% end
% end
% 
% 
% data_path = fileparts(mfilename('fullpath'));
% if datatype==0
%    R_results_path=fullfile(data_path,['Results/synthetic/','n',num2str(n_s),'/R_results']);
%    save(R_results_path,'R_mean','R_variance');
% end
% 
% 
% 
% 




