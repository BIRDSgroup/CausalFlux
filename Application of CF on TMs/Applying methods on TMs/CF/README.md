# Guide for the codes/data in this folder

## CF_MTR 
Contains all the codes/data for applying CF-MTR on the toy models. Please look at the "README.md" file in this folder for further details

## CF_S
Contains all the codes/data for applying CF-S on the toy models. Please look at the "README.md" file in this folder for further details

## Data_generated
Contains the csv files containing the results from applying CF-S, CF-MTR, TRIMER, and GIMME on the toy models. Please look at the "README.md" file in this folder for further details.

The data from this folder is used by the ".R" scripts to generate the figures in the manuscript for the toy models

## Instructions on running the "Sep_act_vs_pred.R" script
Running this script will generate normalized actual vs normalized predicted fluxes plot for the TMs across the WT/KO cases under three different exchange rates. The plots for all the methods (CF-S (FVA and FBA implementation), CF-MTR (FVA and FBA implementation), TRIMER, and GIMME (FVA and FBA implementation)) will be generated.

1) Make sure the working directory in line 115 corresponds to the working diretory where you want the figures to be generated. This will create "Sep_act_vs_pred_plot" folder in that directory you have specified.

2) Make sure the lines 23, 28, and 34 points to the folder "TM1" in the "Data_generated" folder

3) Make sure the lines 153, 158, and 163 points to the folder "TM2" in the "Data_generated" folder

4) Make sure the lines 279, 284, and 289 points to the folder "TM3" in the "Data_generated" folder

## Instructions on running the "SC_sep_plots.R" script
Running this script will generate the SC plots (between the methods), where SC was computed between the actual and predicted fluxes plot for the TMs across the WT/KO cases under three different exchange rates. The plots for all the methods (CF-S (FVA and FBA implementation), CF-MTR (FVA and FBA implementation), TRIMER, and GIMME (FVA and FBA implementation)) will be generated.

1) Make sure the working directory in line 133 corresponds to the working diretory where you want the figures to be generated. This will create "SC_sep_plots" folder in that directory you have specified. 

2) Make sure the lines 36, 44, and 52 points to the folder "TM1" in the "Data_generated" folder

3) Make sure the lines 64, 72, and 80 points to the folder "TM2" in the "Data_generated" folder

4) Make sure the lines 92, 100, and 108 points to the folder "TM3" in the "Data_generated" folder

## Instructions on running the "BM_act_vs_pred_plots.R" script
Running this script will generate the plots between the normalized actual and normalized predicted fluxes for the biomass reactions from the TMs across the WT/KO cases under three different exchange rates taken together. This will generate plots for CF-S (FVA max), CF-MTR (FVA max), TRIMER and GIMME (FVA max).
1) Make sure the working directory in line 87 corresponds to the working diretory where you want the figures to be generated. This will create "BM_ACT_VS_PRED" folder in that directory you have specified.

2) Make sure the working directory in line 12 points to the folder where all the "TM1", "TM2", and "TM3" from the "Data_generated" folder is present


