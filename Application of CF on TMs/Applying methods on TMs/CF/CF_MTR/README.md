# Instructions to run CF-MTR on TM1, 2 and 3

*Make sure all the folders (TM1,2 AND 3), ".m" files, "getMaxTurnOverRate.m", ".mat" files, and the ".R" script  are in the same working directory*

1) Running this script after making all the necessary modifications below will output the steady-state fluxes for TM1, TM2 and TM3 in WT and all gene KO cases under 3.2, 320 and 3200 exchange rates.
2) It will create a folder "RES" containg the necesary KO output folders generated for 3.2, 320 and 3200 case for each TM1, TM2 and TM3.
3) In each of the KO folders for the exchange rates, separate WT and gene folders will get generated.
4) This folder will contain "FVA_to_check.xlsx" and "FBA_to_check.csv" that stores the final iteration values of the steady-state fluxes

## Changes in "CF_MTR_ALL_TMS.R" script
1) Change the working directory in line 3758 to the working directory where all the files and codes from this folder are stored
2) Give your choice of "p" (binariazation percentage) and "maxi" (maximum iterations) in lines 3761 and 3763 

### Changing the working directory in the below mentioned lines pointing to the GRN (structure and paramters learnt)data and metabolic model information 
3) Change the lines in 979, 1018 to working directory pointing to the "grn_3.2" folder inside "TM1" folder
4) Change the lines in 981, 1020 to working directory pointing to the "grn_320" folder inside "TM1" folder
5) Change the lines in 983, 1022 to working directory pointing to the "grn_3200" folder inside "TM1" folder

6) Change the lines in 2214, 2253 to working directory pointing to the "grn_3.2" folder inside "TM2" folder
7) Change the lines in 2216, 2255 to working directory pointing to the "grn_320" folder inside "TM2" folder
8) Change the lines in 2218, 2257 to working directory pointing to the "grn_3200" folder inside "TM2" folder

9) Change the lines in 3480, 3519 to working directory pointing to the "grn_3.2" folder inside "TM3" folder
7) Change the lines in 3482, 3521 to working directory pointing to the "grn_320" folder inside "TM3" folder
8) Change the lines in 3484, 3523 to working directory pointing to the "grn_3200" folder inside "TM3" folder

## Changes in all ".m" scripts (9 ".m" files are there in total)
1) Make sure the "curr_wd" in all these files (except "getMaxTurnOverRate.m") are changed to the working director where the "CF_MTR_ALL_TMS.R" script is present

