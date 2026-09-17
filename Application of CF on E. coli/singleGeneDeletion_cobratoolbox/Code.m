initCobraToolbox(false); % run with initCobraToolbox(false) to avoid update
changeCobraSolver('gurobi', 'all')

%%%%%%%% single-gene KO for E.coli 
%%% 1323 genes
curr_wd = 'D:\work\Integrated_network_model\Ecoli_intg_ntwk\metabolic_aspect\Auto_RUN\Causal_Surgery\Parallel_Runs\CF_S_1';
cd(curr_wd)

curr_wd = 'D:\work\Integrated_network_model\Ecoli_intg_ntwk\metabolic_aspect\Auto_RUN';
cd(curr_wd)
load('iML1515.mat');
ecoli_file = 'iML1515.mat';
eco_mod = readCbModel(ecoli_file);


eco_mod.c(2669)=1;
cd('D:\work\Integrated_network_model\Ecoli_intg_ntwk\metabolic_aspect\Auto_RUN\Causal_Surgery\New_Results\new_data\')
geneTable = readtable('GS_GT_E_genes.csv', 'ReadVariableNames', true);

% Convert to cell array of gene IDs
geneList = geneTable.e_ids;

[grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = ...
    singleGeneDeletion(eco_mod, 'FBA', geneList);



