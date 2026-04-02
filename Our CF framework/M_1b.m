curr_wd = 'D:\work\Integrated_network_model\Git_hub_codes\Ecoli\CF';      
cd(curr_wd)


% some pre-requisites 
initCobraToolbox(false);
changeCobraSolver('gurobi', 'all');


fileName = 'NEW_iML1515_modified_sink.mat';  %%% load the metabolic model 
TM_0 = readCbModel(fileName);

TM_0.c(:) = 0;

cd(curr_wd)
exch_rxns_dt = readtable("Exch_iml1515.csv");   %%% a dataframe with all exchange reactions in the metabolic model and its index

TM_0.lb(exch_rxns_dt.Var1) = 0;


cd(curr_wd)
essen_exch_rxns_dt = readtable("essential_exch_rxns_ids_iml1515.csv");  %%% list of index of the essential exchange reactions in the metabolic model
TM_0.lb(essen_exch_rxns_dt.Var1)= -0.001;

cd(curr_wd)
LBmediaconstraints = readtable("LB_media_constraints_iML1515.csv");  %%% a datframe with index of the exchange reactions for media condition and the upper bounds for those reactions
% LBmediaconstraints = readtable("M9_serine_galacturonate.csv");
% LBmediaconstraints = readtable("TSBmed.csv");


 TM_0.lb(LBmediaconstraints.media_list)=-LBmediaconstraints.media_list_ub;


cd(curr_wd)
TM_0.lb(181) = -readvars("Exch_G.csv");   %%% glucose exchange reaction in the metabolic model
TM_0.lb(1982) = -readvars("Exch_O.csv");  %%% oxygen exchange reaction in the metabolic model


TM_0.ub(2713:2731) = 10; %%% extra sink reactions added based on metabolic feedback information 



[min_x_model_round_0, max_x_model_round_0] = fluxVariability(TM_0);

rxn_abbrev = TM_0.rxns;
minimum_flux_round_0 = num2cell(min_x_model_round_0);
maximum_flux_round_0 = num2cell(max_x_model_round_0);

fva_op_round_0 = [rxn_abbrev,minimum_flux_round_0, maximum_flux_round_0];

cd(curr_wd)
writecell(fva_op_round_0,"FVA_1b_obj_0_P1.xlsx")


TM_1 = TM_0;

TM_1.c(2669)=1; %%% select the biomass reaction for the metabolic model


[min_x_model_round_1, max_x_model_round_1] = fluxVariability(TM_1);

rxn_abbrev = TM_1.rxns;
minimum_flux_round_1 = num2cell(min_x_model_round_1);
maximum_flux_round_1 = num2cell(max_x_model_round_1);

fva_op_round_1 = [rxn_abbrev,minimum_flux_round_1, maximum_flux_round_1];

cd(curr_wd)
writecell(fva_op_round_1,"FVA_to_check_P1.xlsx")


TM_1.c(2669)=1; %%% select the biomass reaction for the metabolic model


sol = optimizeCbModel(TM_1);
solv = num2cell(sol.v);
fba_sol = [rxn_abbrev,solv];
writecell(fba_sol,"FBA_to_check_P1.csv")
