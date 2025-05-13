Hierarchical Bayesian inference for community detection and connectivity of functional brain networks

Version 1.0
12-April-2025
Copyright (c) 2025, Lingbin Bian
---------------------------------------------------------------------------
Main functions of community detection:
CommuDetectLBM.m (individual-level community detection)
CommuDetectGroup.m (group-level community detection)

A simple demonstration:
DEMO_community_detection.m
---------------------------------------------------------------------------
The code of the study in the paper is under the directory: HBI_CommuDetect
The code for processing working memory task-fMRI is under the directory: tfMRI_Data_HCP
Some pictures of the results are under the directory: Pictures_community_labels
---------------------------------------------------------------------------
How to sun the code of the study in the paper: 
---------------------------------------------------------------------------
Experiments: synthetic data analysis
---------------------------------------------------------------------------
1.Use create_simu_dir.sh to create the directory of synthetic data.
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
---------------------------------------------------------------------------
2.Use MANIP_Gaussian_data_generator.m to generate the synthetic data
---------------------------------------------------------------------------
3.Estimate individual community memberships (individual-level analysis)

DEMO_individual_LBM.m (synthetic data: data_type=0)
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

DEMO_individual_modularity.m (synthetic data: data_type=0)
                                 (DIIV)    (SNR)     (gamma)
Results---synthetic_modularity---DIIV10---n0.3162---1 1.2 ... 2
                |                      ---n0.5623
                |                      ---n1
                |                      ---n1.7783
                |                      ---n3.1623
                |                      ---n5.6234
                |                      ---n10
                |                      ---n17.7828
                |----------------DIIV20---n0.3162
                |                      ---n0.5623
                |                      ---n1
                |                      ---n1.7783
                |                      ---n3.1623
                |                      ---n5.6234
                |                      ---n10
                |                      ---n17.7828
                ...


DEMO_individual_multilayer_modularity.m (comparison with multilayer modularity)
---------------------------------------------------------------------------
4.Comparing the estimated individual community mememberships with the ground truth (individual-level analysis) using normalized mutual information (NMI)
MANIP_individual_calcu_NMI_LBM.m 
MANIP_individual_calcu_NMI_modularity.m
MANIP_individual_calcu_NMI_multilayer_modularity.m
---------------------------------------------------------------------------
5.Statistical analysis (individual-level analysis)
DEMO_individual_statistical_analysis_modularity.m (LBM vs modularity)
DEMO_individual_statistical_analysis_multilayer_modularity.m (LBM vs multilayer modularity)
---------------------------------------------------------------------------
6.SNR analysis, NMI of LBM against different levels of SNR (individual-level analysis)
DEMO_individual_statistical_analysis_snr.m
---------------------------------------------------------------------------
7.Estimate group-level community memberships based on hierarchical Bayesian inference (group-level analysis)
DEMO_group_community_detection_HBI.m

SNR analysis (group-level analysis)
MANIP_group_calcu_NMI_LBM.m
DEMO_group_NMI_SNR_LBM.m

DIIV analysis (group-level analysis)
MANIP_group_calcu_NMI_LBM.m
DEMO_group_NMI_DIIV_LBM.m
---------------------------------------------------------------------------
8.Estimate group-level community memberships based on consensus clustering (group-level analysis)
DEMO_group_community_detection_Consensus.m

SNR analysis (group-level analysis)
MANIP_group_calcu_NMI_consensus.m
MANIP_group_calcu_NMI_majorityvote.m
DEMO_group_NMI_SNR_LBM.m
---------------------------------------------------------------------------
9.Estimating group-level mean and variance connectivity

DEMO_group_connectivity_HBI.m

---------------------------------------------------------------------------
Experiments: real data analysis
---------------------------------------------------------------------------
10.Spit-half reproducibility analysis (comparing LBM, modularity, and multilayer modularity)

DEMO_individual_LBM.m (real data: data_type=1)

MANIP_individual_calcu_NMI_real_LBM.m

DEMO_individual_modularity.m (real data: data_type=1)
MANIP_individual_calcu_NMI_real_modularity.m

DEMO_individual_multilayer_modularity.m (real data: data_type=1)
MANIP_individual_calcu_NMI_real_multilayer_modularity.m
---------------------------------------------------------------------------
11.between subject consistency and subject-specific to group consistency

Between subjects:
MANIP_individual_between_subjects_2b0bfix_NMI_real_LBM.m

subject-specific to group:

MANIP_individual_subject2group_2b0bfix_NMI_real_LBM.m

Between conditions:
MANIP_individual_between_conditions_2b0bfix_NMI_real_LBM
---------------------------------------------------------------------------
12.Visualize brain networks by brain netviewer

DEMO_group_brainnet_viewer.m









