% matlab function - Initialization
% part - 1 
curr_wd = 'D:\work\Integrated_network_model\Ecoli_intg_ntwk\metabolic_aspect\Auto_RUN\Causal_Surgery\Parallel_Runs\CF_S_1';
cd(curr_wd)


% some pre-requisites 
% initCobraToolbox(false);
% changeCobraSolver('gurobi', 'all');
load('D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/environment')
restoreEnvironment(environment);

%fileName = 'Ecoli_model_aerobic_sink_no_media.mat';
%fileName = 'Ecoli_sink_no_media_aerobic_iML1515.mat';
fileName = 'NEW_iML1515_modified_sink.mat';
TM_0 = readCbModel(fileName);

% fileName = 'TR_based_modified_ecoli.mat';
% TM_0 = readCbModel(fileName);

TM_0.c(:) = 0;


cd('D:\work\Integrated_network_model\Ecoli_intg_ntwk\metabolic_aspect\Auto_RUN\Causal_Surgery\Exchange_rxns\')
exch_rxns_dt = readtable("Exch_iml1515.csv");

TM_0.lb(exch_rxns_dt.Var1) = 0;


cd('D:\work\Integrated_network_model\Ecoli_intg_ntwk\metabolic_aspect\Auto_RUN\Causal_Surgery\Exchange_rxns\')
essen_exch_rxns_dt = readtable("essential_exch_rxns_ids_iml1515.csv");
TM_0.lb(essen_exch_rxns_dt.Var1)= -0.001;

 cd('D:\work\Integrated_network_model\Ecoli_intg_ntwk\metabolic_aspect\Auto_RUN\Causal_Surgery\CF_S\MM\LB_media\')
 LBmediaconstraints = readtable("LB_media_constraints_iML1515.csv");
% LBmediaconstraints = readtable("M9_serine_galacturonate.csv");
% LBmediaconstraints = readtable("TSBmed.csv");

% 
 TM_0.lb(LBmediaconstraints.media_list)=-LBmediaconstraints.media_list_ub;


 cd(curr_wd)
TM_0.lb(181) = -readvars("Exch_G.csv");   %iML1515
TM_0.lb(1982) = -readvars("Exch_O.csv");


% TM_0.ub(2383:2448) = 1000;

%
TM_0.ub(2713:2731) = 10; % iML1515



[min_x_model_round_0, max_x_model_round_0] = fluxVariability(TM_0);

rxn_abbrev = TM_0.rxns;
minimum_flux_round_0 = num2cell(min_x_model_round_0);
maximum_flux_round_0 = num2cell(max_x_model_round_0);

fva_op_round_0 = [rxn_abbrev,minimum_flux_round_0, maximum_flux_round_0];

cd(curr_wd)
writecell(fva_op_round_0,"FVA_1b_obj_0_P1.xlsx")

% part - 2

TM_1 = TM_0;

% TM_1.c(926)=1;  % for our IAF1260 - Biomass
%  TM_1.c(1005)=1;  % for TR Ecoli - Biomass

% TM_1.c(2267)=1;  % for our IAF1260 - Indole
% TM_1.c(2272)=1;  % for TR Ecoli - Biomass

TM_1.c(2669)=1; % for iML1515


[min_x_model_round_1, max_x_model_round_1] = fluxVariability(TM_1);

rxn_abbrev = TM_1.rxns;
minimum_flux_round_1 = num2cell(min_x_model_round_1);
maximum_flux_round_1 = num2cell(max_x_model_round_1);

fva_op_round_1 = [rxn_abbrev,minimum_flux_round_1, maximum_flux_round_1];

cd(curr_wd)
writecell(fva_op_round_1,"FVA_to_check_P1.xlsx")

% TM_1.c(926)=1;  % for our IAF1260 - Biomass
%  TM_1.c(1005)=1;  % for TR Ecoli - Biomass

% TM_1.c(2267)=1;  % for our IAF1260 - Indole
% TM_1.c(2272)=1;  % for TR Ecoli - Biomass

TM_1.c(2669)=1; % for iML1515

%TM_1.ub(2383:2448) = 1000;


sol = optimizeCbModel(TM_1);
solv = num2cell(sol.v);
fba_sol = [rxn_abbrev,solv];
writecell(fba_sol,"FBA_to_check_P1.csv")




%%% EXTRA STUFFS
% writecell(TM_0.rxns,"Model_RXN_iAF1260.csv")
% writecell(TM_0.rxnNames,"Model_RXN_names_iAF1260.csv")









