# Instructions to run CausalFlux on TMs

*Make sure all the folders, ".m" files and ".mat" files and the ".R" script given here are in the same working directory*

## Changes in "CausalFlux_for_3TMS.R" script
 Change the working directory in line 2 to the working directory where all the files and codes from this folder are stored.

## Changes in all ".m" scripts (9 ".m" files are there in total)
 Make sure the "curr_wd" in all these files are changed to the working directory where all the files and codes from this folder are stored (same as above).

## Key notes
1) Running the "CausalFlux_for_3TMS.R" script after making all the necessary modifications (as mentioned above) will output the neccessary files and plots reported in the manuscript.
2) The folders generated: *KO_data_TM1_3.2*, *KO_data_TM2_3.2*, *KO_data_TM3_3.2*, *KO_data_TM1_320*, *KO_data_TM2_320*, *KO_data_TM3_320*,*KO_data_TM1_3200*, *KO_data_TM2_3200*, *KO_data_TM3_3200*
     2.1) These contain the CausalFlux prediction output files for the gene KOs cases corresponding to the TMs and exchange rate.
     2.2) "FVA_to_check.xlsx" file (for every gene KO case under TMs/exchange rate) contains the final iterations (after convergence) steady-state reaction flux predictions.
3) The folders generated: *SC_sep_plots*, *Sep_act_vs_pred_plot*, *BM_ACT_VS_PRED*
     3.1) Each of these folders contains plots that corresponds to the findings reported in the manuscript.
     3.2) These plots genereted through this run are used to develop Fig. 2, 3 (from main manuscript) and Fig S3, S4 (from supplementary material)
4) The folders generated: *TM1*, *TM2*, *TM3*
     4.1) Each of these folders contains files for the Spearman correlation mentioned in the manuscript
     4.2) These results were used to draw conclusions/results in the manuscript
