# GMPT: A Computational Method to Dissect Colonization Resistance of the Gut Microbiota against Pathogens

## Generate the simulated data in folder 'generate simulated data'
1. The scripts can generate the simulated data as shown in Fig.2 and Supplementary Fig.1.
2. The main script for GLV dynamics to generate the simulated data is `Generate_sim_data_for_GMPT.m`. Please see the comments for the meanings of parameters in the script.
3. Run `GMPT_ALDEx2_for_all.R` and `GMPT_ANCOM_for_all.R` to get DAA results by ALDEx2 and ANCOM, respectively.

## GMPT performance as shown in Fig.2 and Supplementary Fig.1

1. Download the example data from Dropbox: https://www.dropbox.com/sh/cvf46f512byzh1c/AADh_XNUpQsom8xpMjSooaWXa?dl=0
2. The ecological netowrks shown in Fig.2 stored at `network of Fig.2a,b,c.pdf`. 
3. Here we provide the example of synthetic community with N = 30 to perform GMPT and visualize the results in `Fig_2bc_Plot_cohousing_process.m` for Fig.2bc and `Fig2_d2j_GMPT_performance_by_ALDEx2.m` for Fig.2d-j if we used DAA method of ALDEx2.
4. If using DAA method of ANCOM to perform GMPT, please run the script of `FigS1_GMPT_performance_by_ANCOM.m` to get the GMPT results as shown in Supplementary Fig.1.
5. Note that after run `Fig2_d2j_GMPT_performance_by_ALDEx2.m` or `FigS1_GMPT_performance_by_ANCOM.m`, it will generate `N30_ALDEx2_GMPT_EF05.txt`.
