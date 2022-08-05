# GMPT: A Computational Method to Dissect Colonization Resistance of the Gut Microbiota against Pathogens

## Generate the simulated data in folder 'generate simulated data'
1. The scripts can generate the simulated data as shown in Fig.2 and Supplementary Fig.1.
2. The main script for GLV dynamics to generate the simulated data is `Generate_sim_data_for_GMPT.m`. Please see the comments for the meanings of parameters in the script.
3. Run `GMPT_ALDEx2_for_all.R` and `GMPT_ANCOM_for_all.R` to get DAA results by ALDEx2 and ANCOM, respectively.

## GMPT performance as shown in Fig.2 and Supplementary Fig.1

1. Download the example data from Dropbox: https://www.dropbox.com/sh/i0vvkm3kwd2emzw/AACYyU4b4cTnPF6f-rSEuy4ea?dl=0
2. The ecological netowrks shown in Fig.2 stored at `network of Fig.2a,b,c.pdf`. 
3. Here we provide the example of synthetic community with N = 30 to perform GMPT and visualize the results in `Fig_2bc_Plot_cohousing_process.m` for Fig.2bc and `Fig2_d2j_GMPT_performance_by_ALDEx2.m` for Fig.2d-j if we used DAA method of ALDEx2.
4. If using DAA method of ANCOM to perform GMPT, please run the script of `FigS1_GMPT_performance_by_ANCOM.m` to get the GMPT results as shown in Supplementary Fig.1.
5. Note that when running `Fig2_d2j_GMPT_performance_by_ALDEx2.m` or `FigS1_GMPT_performance_by_ANCOM.m`, it will generate `N30_ALDEx2_{***}.txt` and `N30_ANCOM_{***}.txt`. `***` indicates `GMPT`, `CaseControl` (equals `MWAS_1`) and `TwoParts` (equals `MWAS_2`). There files store the 29x3 to perform binary classification of zero or non-zero interaction of other species on pathogen X. Thus, the 1st column is species index (except 1 because it equals pathogen X), 2nd column is effect value of species after performing DAA of ALDEx2 or ANCOM, and 3rd column is ground truth of species interaction on pathogen X.
6. `N30_Binary_classification_metrics_by_ALDEx2.txt` and `N30_Binary_classification_metrics_by_ANCOM.txt` are results of binary classification of zero or non-zero interaction of other species on pathogen X. They are used in `Fig2_d2j_GMPT_performance_by_ALDEx2.m` or `FigS1_GMPT_performance_by_ANCOM.m` to plot Fig.2h,i or Supplementary Fig.1e,f, respectively.
