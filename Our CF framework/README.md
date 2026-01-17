# A walk-through of our CF framework

## Introduction
This file will provide a walk-through of our CF framework. The inputs are same for the CF-S and CF-MTR versions. 

## Calling the function: 
Keep the working directory (curr_dir) to the path containing the matlab scripts and R code.
```R
CF_function <- function(curr_wd,xun,vgval,voval,mi, u, exch_rate, ge,gsl,gpl,mgr,gs,fmr)

```
|Input arguments|Description|
   |---|---|
   |curr_wd|Working directory (directory where ".R", ".m" and ".mat" (metabolic model) scripts are located) |
   |xun|Gene (gene symbol) to be perturbed (can provide one gene or a vector of genes in a for-loop for perturbation)|
   |vgval|Glucose uptake rate (in mmol/gDCW/hr)|
   |voval|Oxygen uptake rate (in mmol/gDCW/hr)|
   |mi|Maximum number of iterations|
   |u|index for gene (1 for one gene or "i" index of the gene being perturbed if for-loop is used|
   |exch_rate|exchange rate for other reactions (default: empty)|
   |ge|Binarized gene expression data (samples x genes)|
   |gsl|GRN structure learning object (bn learn object) |
   |gpl|GRN parameter learning object (bn learn object) |
   |mgr|Feedback from metabolites to genes|
   |gs|Genes and Subsystems|
   |fmr|List of exchange/sink reactions corresponding to the feedback metabolites|



### Formart for specific input arguments
1) **Feedback from metabolites to genes**
   
   The data for this information should of the following format. The columns (column names) must be:
   |Columns|Description|
   |---|---|
   |Metabolites|name of the metabolite that have feedback to the genes|
   |Met_symbol|Symbol of the metabolites in the metabolic model|
   |TF-gene|Gene symbol of TF targeted by metabolites|
   |Interaction|Interaction of metabolite on TF (+ve or -ve regulation)|
   |sink_met_symbol|Sink reactions names|
   
2) **Genes and Subsystems**

    This data should contain:
   |Columns|Description|
   |---|---|
   |Gene| Bigg ids of genes from the metabolic model|
   |Subsystem|Subsystem in which the gene is present|
   |gene_name|Gene symbol |
   
3)  **Metabolic model**

    Metabolic model corresponding to the system being studied. If all metabolites in the **mgr** data are already represented as exchange reactions in the current model, no changes are required. Otherwise, sink reactions should be added for the feedback metabolites.


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



