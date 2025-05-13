function [Timeseries] =load_realts(atlas,subid,session_n,K_seg)
% This function loads the real time series from the subdirectories in Data folder.
% 
% Input: subid: the subject id, eg. '100307'
%        session_n: 1 (session LR), 2 (session RL)
%        K_seg: number of communities
% Output: Timeseries: a struct of time series with parameters
%
% Version 1.0
% 19-Feb-2020
% Copyright (c) 2020, Lingbin Bian
% -------------------------------------------------------------------------
data_path = fileparts(mfilename('fullpath'));
if isempty(data_path), data_path = pwd; end
if atlas==1
   if session_n==1
      subdir_path=fullfile(data_path,'Data/whole_brain_timeseries',subid,'timeseries_WM_LR');
   elseif session_n==2
      subdir_path=fullfile(data_path,'Data/whole_brain_timeseries',subid,'timeseries_WM_RL');
   end
else
   if session_n==1
      subdir_path=fullfile(data_path,'Data/whole_brain_timeseries_Kong',subid,'timeseries_WM_LR');
   elseif session_n==2
      subdir_path=fullfile(data_path,'Data/whole_brain_timeseries_Kong',subid,'timeseries_WM_RL');
   end
end
   load(fullfile(subdir_path,'ROI_whole.txt'));
   signal=ROI_whole';
   [N,T]=size(signal);  % data matrix size
   Timeseries=struct('name','Real',...
                     'signal',signal,...
                     'changepoint_truth',[],...
                     'LatentLabels',[],...
                     'T',T,...
                     'N',N,...
                     'K_seg',K_seg);

end