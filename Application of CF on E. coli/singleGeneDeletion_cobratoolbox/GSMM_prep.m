initCobraToolbox(false); % run with initCobraToolbox(false) to avoid update
changeCobraSolver('gurobi', 'all')

cd("/work/Integrated_network_model/BS/Req_mats/GSMM/")
load('iYO844.mat');

bs_mod = iYO844;


extra_mets = {'fdp_c',...
             'mqn7_c',...
             'fe2_c',...
             'h2o2_c',...
             'o2_c',...
             'mg2_c',...
             'zn2_c',...
             'gln__L_c',...
             'ile__L_c',...
             'val__L_c',...
             'alltn_c',...
             'akg_c',...
             'ac_c',...
             'alltt_c',...
             'cit_c',...
             'galur_c',...
             'sbt__D_c',...
             'xyl__D_c',...
             'arab__L_c',...
             'arg__L_c',...
             'citr__L_c',...
             'glu__L_c',...
             'his__L_c',...
             'orn__L_c',...
             'urate_c',...
             'tre6p_c',...
             '12dhg3p_c',...
             'malcoa_c',...
             'acser_c',...
};



ex_rxn = find(findExcRxns(bs_mod));
ex_rnx_names  = bs_mod.rxnNames(ex_rxn);
ex_rnx_sym  = bs_mod.rxns(ex_rxn);

% Convert both to string arrays for convenience
A = extra_mets;     % 1x66
B = ex_rnx_sym;     % 337x1

match_indices = zeros(size(A));

% Loop over elements of A
for i = 1:length(A)
    % Get base name (before underscore)
    core_name = regexp(A{i}, '^[^_]+', 'match', 'once');
    
    % Find matching index in B (case-sensitive)
    match_found = false;
    for j = 1:length(B)
        if contains(B{j}, core_name)
            match_indices(i) = j;  % store first match
            match_found = true;
            break;
        end
    end
    
    % If no match found, keep it as 0 (already initialized)
end

t_m = transpose(match_indices);

extracted_values = cell(size(t_m));

% Loop through t_m
for i = 1:length(t_m)
    idx = t_m(i);
    if idx >= 1 && idx <= length(ex_rnx_sym) && mod(idx,1) == 0
        extracted_values{i} = ex_rnx_sym{idx};  % valid index
    else
        extracted_values{i} = 0;       % for idx == 0 or invalid
    end
end

t_em = transpose(extra_mets);

%% modified extra mets 

extra_mets_mets = {'fdp_c',...
             'mqn7_c',...
             'alltt_c',...
             'tre6p_c',...
             '12dhg3p_c',...
             'malcoa_c',...
             'acser_c',...
};

bs_mod_mod = bs_mod;
bs_mod_mod.c(:)=0;

bs_mod_mod = addSinkReactions(bs_mod_mod,extra_mets_mets,[zeros(numel(extra_mets_mets),1)],[zeros(numel(extra_mets_mets),1)+1000]);

cd("/work/Integrated_network_model/BS/Req_mats/GSMM/")
save("iYO844_new_model.mat","bs_mod_mod")


%%%
%%%
bs_mod.c(150) = 1;
[grRatio, grRateKO, grRateWT, rxnList] = singleRxnDeletion(bs_mod);
essentialRxns = rxnList(grRatio == 0);
essentialRxns_names = bs_mod.rxns(rxnList);
%%%


% Single-gene KO in BS 
% Read gene list

curr_wd = 'D:\work\Integrated_network_model\BS\parallel_runs\CF_S_1';
cd(curr_wd)
fileName = 'iYO844_new_model.mat';
TM_0 = readCbModel(fileName);

curr_wd = 'D:\work\Integrated_network_model\BS\Req_mats\KO_info\';
cd(curr_wd)

TM_0.c(332)=1;
%geneTable = readtable('SKO_BS.csv', 'ReadVariableNames', false);
geneTable = readtable('GSMM_gene_KO_BS.csv', 'ReadVariableNames', false);

% Convert to cell array of gene IDs
geneList = geneTable.Var1;



cd('D:\work\Integrated_network_model\BS\Req_mats\GSMM\')
exch_rxns_dt = readtable("iYO844_all_exchnage_rxns.csv");

TM_0.lb(exch_rxns_dt.Var1) = 0;


cd('D:\work\Integrated_network_model\BS\Req_mats\GSMM\')
essen_exch_rxns_dt = readtable("iYO844_essential_exch_rxns_ids.csv");
TM_0.lb(essen_exch_rxns_dt.Var1)= -0.001;

cd('D:\work\Integrated_network_model\BS\Req_mats\GSMM\')
LBmediaconstraints = readtable("LBmed_iYO844.csv");

 TM_0.lb(LBmediaconstraints.media_list)=-LBmediaconstraints.media_list_ub;


 cd('D:\work\Integrated_network_model\BS\parallel_runs\CF_S_1')
TM_0.lb(59) = -readvars("Exch_G.csv");   %BS
TM_0.lb(359) = -readvars("Exch_O.csv");

%
TM_0.ub(1251:1257) = 10; % BS



% Perform single gene deletions
[grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = ...
    singleGeneDeletion(TM_0, 'FBA', geneList);



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



