#!bin/bash
mkdir ../whole_brain__timeseries
cd ../subjects
for subj in `cat subject.txt` ; do

   mkdir ../whole_brain_timeseries/$subj
   mkdir ../whole_brain_timeseries/$subj/timeseries_WM_LR
   mkdir ../whole_brain_timeseries/$subj/timeseries_WM_RL
done
cd ../Shell_script_control

cd ../subjects
for subj in `cat subject.txt` ; do
   fslmeants -i $subj/MNINonLinear/Results/tfMRI_WM_LR/tfMRI_WM_LR.nii.gz --label=../Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm.nii.gz -o ../whole_brain_timeseries/$subj/timeseries_WM_LR/ROI_whole.txt  
   fslmeants -i $subj/MNINonLinear/Results/tfMRI_WM_RL/tfMRI_WM_RL.nii.gz --label=../Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm.nii.gz -o ../whole_brain_timeseries/$subj/timeseries_WM_RL/ROI_whole.txt  
   
done 
cd ../Shell_script_control

