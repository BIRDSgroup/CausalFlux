# Instructions for Ablation analysis 

## To run the ablated model 

### Option A: To run CausalFlux from scratch for the ablation models (5 ablation models; each perform 798 single-gene KOs)
1) Use the code from "CausalFlux_Runs" section [Go to the CausalFlux_Runs](../CausalFlux_Runs)
2) Make the modification related to ablation as explained below
   
    2.1) For "crp" ablation: Add the below lines of code in line 1180 of "CausalFlux.R" script
   ```R
   not_req_MGR_genes <- c("crp")
   xi <- which(MGR$`TF-gene` %in% not_req_MGR_genes)
   MGR <- MGR[-xi,]
   ```

   2.2) For "fur" ablation: Add the below lines of code in line 1180 of "CausalFlux.R" script
   ```R
   not_req_MGR_genes <- c("fur")
   xi <- which(MGR$`TF-gene` %in% not_req_MGR_genes)
   MGR <- MGR[-xi,]
   ```
   2.3) For "cra" ablation: Add the below lines of code in line 1180 of "CausalFlux.R" script
   ```R
   not_req_MGR_genes <- c("cra")
   xi <- which(MGR$`TF-gene` %in% not_req_MGR_genes)
   MGR <- MGR[-xi,]
   ```
   2.4) For "10 metabolic feedback genes" ablation: Add the below lines of code in line 1180 of "CausalFlux.R" script
   ```R
   not_req_MGR_genes <- c("crp" , "fur" , "cra" , "lrp",  "pdhR", "argR" ,"purR" ,"cysB", "fhlA", "nagC")
   xi <- which(MGR$`TF-gene` %in% not_req_MGR_genes)
   MGR <- MGR[-xi,]
   ```
   2.5) For "10 metabolic feedback genes" ablation: Add the below lines of code in line 1180 of "CausalFlux.R" script
   ```R
   not_req_MGR_genes <- c("crp" , "fur" , "cra" , "lrp",  "pdhR", "argR" ,"purR" ,"cysB", "fhlA", "nagC",
                       "metJ" ,"paaX" ,"cytR", "argP", "cbl" , "glpR", "tyrR","malT", "gntR", "araC",
                       "trpR", "galR", "galS" ,"allR","exuR" ,"deoR", "mhpR" ,"lsrR" ,"hcaR", "idnR" ,
                       "uxuR" ,"rbsR", "xylR", "fucR","prpR" ,"caiF", "metR" ,"rhaS" ,"nanR" ,"rutR")
   xi <- which(MGR$`TF-gene` %in% not_req_MGR_genes)
   MGR <- MGR[-xi,]
   ```
   
### Option B: Use the pre-computed predictions from the ablation models
1) "Ecoli_Ablation_Predictions.csv" already contains the predictions from the CausalFlux model run on the 5 different ablation cases for 798 single-gene KOs
2) The given drive link contains the folders/files with the CausalFlux predictions on the 5 different ablation cases for 798 single-gene KOs. "Ecoli_Ablation_Predictions.csv" was generated using this data.
   
## Use "Ablation.R" to plot the results
Change the "setwd()" in line 5 to match the working directory where the "Ablation.R", "Ecoli_Ablation_Predictions.csv", "Metabolite_Gene_Regulation.RDS" and "bn_TR.RDS" are stored 

## Outputs
Running the script generates plots and CSV files used for Figure 5, Supplementary Figures S5, Supplementary Tables 4 and 5.


