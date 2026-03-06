Hierarchical Bayesian inference for community detection and connectivity of functional brain networks
---------------------------------------------------------------------------

Version 1.0

6-March-2026

Copyright (c) 2026, Lingbin Bian

Main functions of community detection:
---------------------------------------------------------------------------

CommuDetectLBM.m

CommuDetectGroup.m

A simple demonstration:
---------------------------------------------------------------------------

DEMO_community_detection.m

The code of the study in the paper: 
---------------------------------------------------------------------------

HBI_CommuDetect

The code for processing working memory task-fMRI:
---------------------------------------------------------------------------

tfMRI_Data_HCP
---------------------------------------------------------------------------

How to run the code of the study in the paper: 
---------------------------------------------------------------------------

Experiments: synthetic data analysis
---------------------------------------------------------------------------

1.Use create_simu_dir.sh to create the directory of synthetic data.
---------------------------------------------------------------------------

Directory structure:

                   (DIIV)    (SNR)     (subject)
Data---synthetic---DIIV10---n0.3162---1001...1100
                |        ---n0.5623
                |        ---n1
                |        ---n1.7783
                |        ---n3.1623
                |        ---n5.6234
                |        ---n10
                |        ---n17.7828
                |--DIIV20---n0.3162
                |        ---n0.5623
                |        ---n1
                |        ---n1.7783
                |        ---n3.1623
                |        ---n5.6234
                |        ---n10
                |        ---n17.7828
                ...

where the folder name 'DIIV' indicates the degree of inter-individual variation, 
'n' indicates the standard deviation of the noise.


2.Use MANIP_Gaussian_data_generator.m to generate the synthetic data
---------------------------------------------------------------------------

3.Estimate individual community memberships (individual-level analysis)
---------------------------------------------------------------------------

ANAL_individual_LBM.m (synthetic data: data_type=0)

Directory structure:

                          (DIIV)    (SNR) 
Results---synthetic_LBM---DIIV10---n0.3162
                |               ---n0.5623
                |               ---n1
                |               ---n1.7783
                |               ---n3.1623
                |               ---n5.6234
                |               ---n10
                |               ---n17.7828
                |---------DIIV20---n0.3162
                |               ---n0.5623
                |               ---n1
                |               ---n1.7783
                |               ---n3.1623
                |               ---n5.6234
                |               ---n10
                |               ---n17.7828
                ...


4.Comparing LBM, (multilayer) modularity (individual-level analysis)
---------------------------------------------------------------------------

ANAL_Bayesian_optimization_LBM.m

ANAL_Bayesian_optimization_modularity.m

ANAL_Bayesian_optimization_multilayer_modularity.m

ANAL_individual_LBM.m

ANAL_individual_modularity.m

ANAL_individual_multilayer_modularity.m

ANAL_LBM_vs_modularity_opt.m

5.SNR analysis, NMI of LBM against different levels of SNR (individual-level analysis)
---------------------------------------------------------------------------

ANAL_individual_statistical_analysis_snr.m

6.Estimate group-level community memberships based on hierarchical Bayesian inference (group-level analysis)
---------------------------------------------------------------------------

ANAL_group_community_detection_HBI.m

SNR analysis (group-level analysis)

ANAL_group_calcu_NMI_LBM.m

ANAL_group_NMI_SNR_LBM.m

DIIV analysis (group-level analysis)

ANAL_group_calcu_NMI_LBM.m

ANAL_group_NMI_DIIV_LBM.m


7.Estimate group-level community memberships
---------------------------------------------------------------------------

SNR analysis (group-level analysis)

ANAL_group_calcu_NMI_LBM.m

ANAL_group_calcu_NMI_majorityvote.m

ANAL_group_NMI_SNR_LBM.m

8.Estimating group-level mean and variance connectivity
---------------------------------------------------------------------------
ANAL_group_connectivity_HBI.m


Experiments: real data analysis
---------------------------------------------------------------------------

9.Spit-half reproducibility (comparing LBM, modularity, and multilayer modularity)
---------------------------------------------------------------------------

ANAL_Bayesian_optimization_LBM_real.m

ANAL_Bayesian_optimization_modularity_real.m

ANAL_Bayesian_optimization_multilayer_modularity_real.m

ANAL_individual_LBM_real.m

ANAL_individual_modularity_real.m

ANAL_individual_multilayer_modularity_real.m


10.between subject consistency and subject-specific to group consistency
---------------------------------------------------------------------------

Between subjects:

ANAL_individual_between_subjects_2b0bfix_NMI_real_LBM.m

subject-specific to group:

ANAL_individual_subject2group_2b0bfix_NMI_real_LBM.m

Between conditions:

ANAL_individual_between_conditions_2b0bfix_NMI_real_LBM


11.Visualize brain networks by brain netviewer
---------------------------------------------------------------------------

ANAL_group_brainnet_viewer.m








