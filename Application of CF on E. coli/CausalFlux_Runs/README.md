# Instuctions to run single-gene KO or multiple single-gene KOs in *E. coli* 

*Make sure all the folders, ".m" files, ".csv" files, ".mat" file, and the ".R" script given here are in the same working directory*

## Changes to be made in "CausalFlux.R" script
1) Change the working directory argument "curr_wd" in the line 1160. Make sure the  working directory where all the files and codes from this folder are stored.
2) Example for running the "CausalFlux_runs" function in two cases are given here (line 1190 and 1195). Choose any one and comment out the other for your objective
3) Modify the "ips" argument accordingly
   
## Changes to be made in three ".m" scripts 
1) Make sure the working directory "curr_wd" (line 1) in all these scripts are changed to the working directory where all the files and codes from this folder are stored (same as above).
2) By default (line 25) "LB media" is the environment for the single-gene KOs. You can change this other media conditions like "Minimal M9 media" (line 26) or "TSB media" (line 27) depending on your choice

## Key Notes 

1) Running the "CausalFlux.R" script after making all the necessary modifications (as mentioned above) will output the neccessary results which have been reported in the manuscript
2) You can run this script under two cases
   
2.1) Case - 1: This setting is to run our CausalFlux algorithm for only one single-gene KO or WT condition

a) Make sure the "ips" argument is vector of length 1.

b) For example: ips <- c("serC") or ips <- c("WT") or ips <- c("gabT")  

2.2) Case - 2: This setting is to run our CausalFlux algorithm for multiple single-gene KOs (including WT condition as well)

a) Make sure the "ips" argument is vector of length at least 2.

b) For example: ips <- c("WT,"serC") or ips <- c("serC", "gabT", "argD")

c) After running, the folder **CausalFlux_multi_runs** will be generated containing the prediction results for each gene KO in a separate folder

3) The excel file *FVA_to_check_P1.xlsx* under the folders corresponding to gene KO contains the results after the final iteration/convergenece of algorithm. The biomass reaction "BIOMASS_Ec_iML1515_core_75p37M" from this file is looked at for the final result.
4) The default values for the arguments like **maxiteration**, **Glucose_exchange**, **Oxygen_exchange** and **extra_rxn_exchange** are already set here

## 798 single-gene KOs
The predictions from CausalFlux for these 798 single-gene KOs are given [here](https://drive.google.com/drive/folders/1fqiBEupVbZp4vPJ7xnCCexDrCc1OJaQM?usp=sharing). Each of these contains separate files that store different outputs of the algorithm. "FVA_to_check_P1.xlsx" in each of these folders contain the final iteration fluxes values.

The [Overall_Ecoli_preds.csv](https://github.com/BIRDSgroup/CausalFlux/blob/main/Application%20of%20CF%20on%20E.%20coli/CausalFlux_Comparison/Overall_Ecoli_preds.csv) file contains the CausalFlux predictions on the 798 single-gene KOs which has been derived from the individual runs present in the above mentioned Google drive.

## To run for 798 single-gene KOs (faster way)
 1) The list of 798 gene names can be obtained from
 2) For faster runs of CausalFlux you can create multiple folders (K) of "CausalFlux_Runs" with its contents
 3) Open K terminals and run the "CausalFlux.R" scripts from each of these K folders
 4) Each of these K executions of "CausalFlux.R" can be run of "case 2" setting to run for multiple single-gene KOs

       4a) For example: Split 798 single-genes into 12 (K) parts (roughly ~67 KOs)

       4b) Create 12 (K) folders of "CausalFlux_Runs". Open 12 (K) terminals (set working directory accordingly). Run "CausalFlux.R" from these 12 (K) terminals

       4c) Peform 67 single-gene KOs using the "case 2" setting of the "CausalFlux_runs" function in each of these scripts. That is the "ips" would be a 67 length vector  in each of these scripts. 







