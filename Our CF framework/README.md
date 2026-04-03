# A walk-through of our CasusalFlux framework

## Introduction
This file will provide a walk-through of our CasusalFlux framework

## Calling the function: 
Keep the working directory (curr_dir) to the path containing the matlab scripts and R code.
```R
CausalFlux_runs <- function(case,curr_wd,ips,maxiteration,Glucose_exchange,Oxygen_exchange,extra_rxn_exchange, GE, GSL, GPL,MGR,GS,FMR)
```
|Input arguments|Description|
   |---|---|
   |case|1 for scenario to run single-gene KO or WT condition in E. coli; 2 for scenario to run multiple single-gene KOs in E. coli |
   |curr_wd|Working directory (directory where ".R", ".m" and ".mat" (metabolic model) scripts and other necessary files/folders are located)|
   |ips|Gene (gene symbol) to be perturbed (can provide one gene or a vector of genes in a for-loop for perturbation)|
   |maxiteration|Maximum number of iterations|
   |Glucose_exchange|Glucose uptake rate (in mmol/gDCW/hr)|
   |Oxygen_exchange|Oxygen uptake rate (in mmol/gDCW/hr)|
   |extra_rxn_exchange|exchange rate for other reactions (default: empty)|
   |GE|Binarized gene expression data (samples x genes)|
   |GSL|GRN structure learning object (bn learn object) |
   |GPL|GRN parameter learning object (bn learn object) |
   |MGR|Feedback from metabolites to genes|
   |GS|Genes and Subsystems|
   |FMR|List of exchange/sink reactions corresponding to the feedback metabolites|



### Format for specific input arguments
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

    Metabolic model corresponding to the system being studied. If all metabolites in the **MGR** data are already represented as exchange reactions in the current model, no changes are required. Otherwise, sink reactions should be added for the feedback metabolites using the "addSinkReactions()" function from COBRA Toolbox.
    An example case is shown below:
```matlab
metabolites_to_be_added_as_sink = {'metA',...
'metB',...
'metC',...
};

metabolic_model_with_sink = addSinkReactions(model = metabolic_model, metabolites = metabolites_to_be_added_as_sink, lb = [zeros(numel(metabolites_to_be_added_as_sink),1)], ub = [zeros(numel(metabolites_to_be_added_as_sink),1)+1000])
```

### (B) To be changed in the ".m" scripts

#### M_1b.m

|Line|Description|Notes|
|---|---|---|
|Line 1|Working directory|directory where ".R", ".m" and ".mat" (metabolic model) scripts and other necessary files/folders are located|
|Line 9|Metabolic model|Metabolic model with added sink reactions|
|Line 13|Exchnage reactions of metabolic model| a dataframe with index of all exchange reactions in the metabolic model and its name|
|Line 19|Essential exchnage reactions of metabolic model| a vector with index of only essential exchange reactions in the metabolic model|
|Line 23|Media constraints| a dataframe with index of the exchange reactions for media condition and the upper bounds for those reactions|
|Line 31|Glucose reaction|index of glucose exchange reaction|
|Line 32|Oxygen reaction|index of oxygen exchange reaction|
|Line 34|sink reaction indices|indices specifying start and end of added sink reactions|
|Line 52|Biomass reaction|index of biomass reaction|
|Line 67|Biomass reaction|index of biomass reaction|

#### M_2a.m

|Line|Description|Notes|
|---|---|---|
|Line 1|Working directory|directory where ".R", ".m" and ".mat" (metabolic model) scripts and other necessary files/folders are located|
|Line 9|Metabolic model|Metabolic model with added sink reactions|
|Line 13|Exchnage reactions of metabolic model| a dataframe with index of all exchange reactions in the metabolic model and its name|
|Line 19|Essential exchnage reactions of metabolic model| a vector with index of only essential exchange reactions in the metabolic model|
|Line 23|Media constraints| a dataframe with index of the exchange reactions for media condition and the upper bounds for those reactions|
|Line 31|Glucose reaction|index of glucose exchange reaction|
|Line 32|Oxygen reaction|index of oxygen exchange reaction|
|Line 34|sink reaction indices|indices specifying start and end of added sink reactions|
|Line 38|Biomass reaction|index of biomass reaction|
|Line 58|Biomass reaction|index of biomass reaction|

#### M_2d.m

|Line|Description|Notes|
|---|---|---|
|Line 1|Working directory|directory where ".R", ".m" and ".mat" (metabolic model) scripts and other necessary files/folders are located|
|Line 9|Metabolic model|Metabolic model with added sink reactions|
|Line 13|Exchnage reactions of metabolic model| a dataframe with index of all exchange reactions in the metabolic model and its name|
|Line 19|Essential exchnage reactions of metabolic model| a vector with index of only essential exchange reactions in the metabolic model|
|Line 23|Media constraints| a dataframe with index of the exchange reactions for media condition and the upper bounds for those reactions|
|Line 31|Glucose reaction|index of glucose exchange reaction|
|Line 32|Oxygen reaction|index of oxygen exchange reaction|
|Line 34|sink reaction indices|indices specifying start and end of added sink reactions|




