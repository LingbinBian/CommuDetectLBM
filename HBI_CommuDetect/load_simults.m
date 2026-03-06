function [Timeseries,latent_seg,latent_seg_sub,K_seg] =load_simults(simul_id,n_s,vari,hrf_ind)
% This function loads the simulated time series.
%
% Input: simul_id: simulated data id. eg.'1001'
% Output: Timeseries: A struct of time series with parameters

% Version 1.0 
% Copyright (c) 2020, Lingbin Bian
% 19-Feb-2020
% -------------------------------------------------------------------------
data_path = fileparts(mfilename('fullpath'));
if isempty(data_path), data_path = pwd; end
% -------------------------------------------------------
  
subdir_path=fullfile(data_path,['Data/synthetic/DIIV',num2str(vari),'/n',num2str(n_s)],simul_id);
      
load(fullfile(subdir_path,'timeseries.mat'));    

% A struct 'Timeseries' is contained in the workspace
end
