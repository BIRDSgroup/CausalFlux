# Instructions to run CausalFlux on TMs

*Make sure all the folders, ".m" files and ".mat" files and the ".R" script given here are in the same working directory*

## Changes in "CausalFlux_for_3TMS.R" script
 Change the working directory in line 2 to the working directory where all the files and codes from this folder are stored.

## Changes in all ".m" scripts (9 ".m" files are there in total)
 Make sure the "curr_wd" in all these files are changed to the working directory where all the files and codes from this folder are stored (same as above).

## Key notes
1) Running the "CausalFlux_for_3TMS.R" script after making all the necessary modifications (as mentioned above) will output the necessary files and plots reported in the manuscript.
2) The folders generated: *KO_data_TM1_3.2*, *KO_data_TM2_3.2*, *KO_data_TM3_3.2*, *KO_data_TM1_320*, *KO_data_TM2_320*, *KO_data_TM3_320*,*KO_data_TM1_3200*, *KO_data_TM2_3200*, *KO_data_TM3_3200*
   
     2.1) These contain the CausalFlux prediction output files for the gene KOs cases corresponding to the TMs and exchange rate.
   
     2.2) "FVA_to_check.xlsx" file (for every gene KO case under TMs/exchange rate) contains the final iterations (after convergence) steady-state reaction flux predictions.
   
3) The folders generated: *SC_sep_plots*, *Sep_act_vs_pred_plot*, *BM_ACT_VS_PRED*
   
     3.1) Each of these folders contains plots that correspond to the findings reported in the manuscript.

      |Folders|Description|
   |---|---|
   |SC_sep_plots|Plots from Fig 2D, Fig 3[A-D], Fig S3|
   |BM_ACT_VS_PRED|Plots from Fig 3[E-G]|
   |Sep_act_vs_pred_plot|Plots from Fig 2[B-C]; such plots for all the 39 TM cases can be found here|

   "rho_iter1_iter2.pdf" in the working directory is Fig S4
   
     3.2) These plots generated through this run are used to develop Fig. 2, 3 (from the main manuscript) and Fig S3 and S4 (from the supplementary material)
   
3) The folders generated: *TM1*, *TM2*, *TM3*
     |Folders|Description|File names|
   |---|---|---|
   |TM1|Predictions of Actual, CausalFlux (FVAm min, FVA max and FBA), GIMME (FVA min, FVA max and FBA); Spearman correlation between these methods and actual - for all TM1 cases|Actual_Pred_data_3.2_TM_1, Actual_Pred_data_320_TM_1, Actual_Pred_data_3200_TM_1; Cor_PV_RMSE_data_3.2_TM_1, Cor_PV_RMSE_data_320_TM_1, Cor_PV_RMSE_data_3200_TM_1  |
   |TM2|Predictions of Actual, CausalFlux (FVAm min, FVA max and FBA), GIMME (FVA min, FVA max and FBA); Spearman correlation between these methods and actual - for all TM2 cases|Actual_Pred_data_3.2_TM_2, Actual_Pred_data_320_TM_2, Actual_Pred_data_3200_TM_2; Cor_PV_RMSE_data_3.2_TM_2, Cor_PV_RMSE_data_320_TM_2, Cor_PV_RMSE_data_3200_TM_2 |
   |TM3|Predictions of Actual, CausalFlux (FVAm min, FVA max and FBA), GIMME (FVA min, FVA max and FBA); Spearman correlation between these methods and actual - for all TM3 cases|Actual_Pred_data_3.2_TM_3, Actual_Pred_data_320_TM_3, Actual_Pred_data_3200_TM_3; Cor_PV_RMSE_data_3.2_TM_3, Cor_PV_RMSE_data_320_TM_3, Cor_PV_RMSE_data_3200_TM_3 |

     4.1) Each of these folders contains files for the Spearman correlation mentioned in the manuscript
   
     4.2) These results were used to draw conclusions/results in the manuscript

     4.3) Supplementary Data 1: all "Cor_PV_RMSE" data combined from TM1, TM2 and TM3 folders

   ## Data folders
   1) GRN folder contains the learnt GRN structure and parameter information (".rds" files), gene expression data (".csv" file), MGR data, GS data  for each TM and exchange rate
   2) ODE_simu folder contains the simulated fluxes (actual) for TM, WT/KO condition, and exchange rate
