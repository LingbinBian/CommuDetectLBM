#!bin/bash
mkdir ../whole_brain_timeseries_Kong
cd ../subjects
for subj in `cat subject.txt` ; do

   mkdir ../whole_brain_timeseries_Kong/$subj
   mkdir ../whole_brain_timeseries_Kong/$subj/timeseries_WM_LR
   mkdir ../whole_brain_timeseries_Kong/$subj/timeseries_WM_RL
done
cd ../Shell_script_control

cd ../subjects
for subj in `cat subject.txt` ; do
   fslmeants -i $subj/MNINonLinear/Results/tfMRI_WM_LR/tfMRI_WM_LR.nii.gz --label=../200Parcels_Kong2022_17Networks_FSLMNI152_2mm.nii.gz -o ../whole_brain_timeseries_Kong/$subj/timeseries_WM_LR/ROI_whole.txt  
   fslmeants -i $subj/MNINonLinear/Results/tfMRI_WM_RL/tfMRI_WM_RL.nii.gz --label=../200Parcels_Kong2022_17Networks_FSLMNI152_2mm.nii.gz -o ../whole_brain_timeseries_Kong/$subj/timeseries_WM_RL/ROI_whole.txt  
   
done 
cd ../Shell_script_control

