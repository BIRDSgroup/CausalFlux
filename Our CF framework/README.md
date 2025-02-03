# A walk-through of our CF framework

## Introduction
This file will provide a walk-through of our CF framework. The inputs are same for the CF-S and CF-MTR versions. 

## Inputs for the CF framework: 
### (A) To be changed in the ".R" scripts
1) **Reconstructed GRN and the parameters learned from this GRN**
   
   GRN reconstruction and parameter leazrning can be done independently, and should be stored in the folder **GRN_REQ_WO_BIGG**. The working directory of this folder must be set in the ".R" scripts for CF-S and CF-MTR.
   The "bn" objects corresponding to these data can be stored as RDS data.
2) **Binarized gene expression data**
   
   The binarized gene expression data used to reconstruct/learn the parameters of the GRN. This should also be in the folder **GRN_REQ_WO_BIGG**. The working directory of this folder must be set in the ".R" scripts for CF-S 
   and CF-MTR
3) **Feedback from metabolites to genes**
   
   The data for this information should of the following format. The columns (column names) must be:
   |Columns|Description|
   |---|---|
   |Metabolites|name of the metabolite that have feedback to the genes|
   |Met_symbol|Symbol of the metabolites in the metabolic model|
   |TF-gene|Gene symbol of TF targeted by metabolites|
   |Interaction|Interaction of metabolite on TF (+ve or -ve regulation)|
   |sink_met_symbol|Sink reactions names|

   This data (as RDS data) should be stored in the folder **MM_REQ**. The working directory of this folder must be set in the ".R" scripts for CF-S and CF-MTR.
   
4) **Genes and Subsystems**

    This data should contain:
   |Columns|Description|
   |---|---|
   |Gene| Bigg ids of genes from the metabolic model|
   |Subsystem|Subsystem in which the gene is present|
   |gene_name|Gene symbol |

   This data (as RDS data) should be stored in the folder **MM_REQ**. The working directory of this folder must be set in the ".R" scripts for CF-S and CF-MTR.
   
### (B) To be changed in the ".m" scripts
1)  **Metabolic model**

    Metabolic model corresponding to the system being studied. In CF-S, add sink reactions (extra reactions added to metabolic model) for the metabolites that have feedback to the genes in GRN. In the case of CF-MTR, keep the metabolic model as such (without sink reactions). Make sure to load the correct metabolic models in all the three ".m" scripts.

### (C) To be given as arguments to the CF function

This "CF_function" (see below) is called as "CF_S_Func" in the case of CF-S and "CF_MTR_Func" in the case of CF-MTR. Nevertheless, the inputs to these functions are still the same.
```R
CF_function(curr_wd, xun, pe, vgval, voval, mi, u)
```
Inputs: 

(a) curr_wd: working directory where all the ".R", ".m", ".mat" (metabolic model) and other folders like "GRN_REQ_WO_BIGG", "MM_REQ" are present. This is also where the results will be stored

(b) xun: Gene (gene symbol) to be perturbed/knocked out

(c) pe: Binarization percentage

(d) vgval: Glucose uptake value (in mmol/gDCW/hr)

(e) voval: Oxygen uptake value (in mmol/gDCW/hr)

(f) mi: Maximum number of iterations

(g) u: It corresponds to the index of the gene being perturbed. 


## Instructions with respect to CF-S:
### (A) Need to changed in the ".R" script

1) Change lines 1148, 1155 to reflect the changes for A.1 and A.2 above
   
2) Change lines 1176 to reflect the changes for A.3 and A.4 above

3) Change line 1318 to set the start and end indices of the sink reactions added to the model. This is input to the "to_check_sink_rxn_iter" function.

### (B) To be changed in the ".m" scripts

#### M_1b.m

1) Change line 3 for current working directory where the codes/data for running the whole CF framework is stored
   
2) Change line 13 for metabolic model name

3) Change line 22 and 23 for changing the indices of glucose and oxygen exchange reactions in the metabolic model (will vary in different metabolic models)

4) Change line 25 to modify the start and end indices of the sink reactions

5) Change line 42, 56 to index number corresponding to the biomass reaction in the metabolic model

#### M_2a.m

1) Change line 2 for current working directory where the codes/data for running the whole CF framework is stored

2) Change line 15 for metabolic model name

3) Change line 25 and 26 for changing the indices of glucose and oxygen exchange reactions in the metabolic model (will vary in different metabolic models)

4) Change line 28 to modify the start and end indices of the sink reactions

5) Change line 32, 51 to index number corresponding to the biomass reaction in the metabolic model

#### M_2d.m

1) Change line 8 for current working directory where the codes/data for running the whole CF framework is stored

2)  Change line 15 for metabolic model name

3)  Change line 24 and 25 for changing the indices of glucose and oxygen exchange reactions in the metabolic model (will vary in different metabolic models)

4)  Change line 27 to modify the start and end indices of the sink reactions

### (C) To be given as arguments to the CF function (in the ".R" script)

1) Check line 1386 when the function is called for performing single gene KO
   
2) Check line 1390 when the function is called for performing multiple single gene KOs. Each KO will be done one by one.




