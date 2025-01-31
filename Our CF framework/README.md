# A walk-thorugh of our CF framework

## Introduction
This file will provide a walk-through of our CF framework. The inputs are same for the CF-S and CF-MTR versions. 

## Inputs for the CF framework: 
### (A) Need to changed in the ".R" scripts
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
5)  **Metabolic model**

    Metabolic model corresponding to the system being studied. In CF-S, add sink reactions (extra reactions added to metabolic model) for the metabolites that have feedback to the genes in GRN. In the case of CF-MTR, keep the metabolic model as such (without sink reactions)
### (B) Need to be given as arguments to the CF function

1) 
