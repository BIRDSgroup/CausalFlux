% matlab function - Initialization
% part - 1 
curr_wd = 'D:\work\Integrated_network_model\Ecoli_intg_ntwk\metabolic_aspect\Auto_RUN\Causal_Surgery\Parallel_Runs\CF_S';
cd(curr_wd)


% some pre-requisites 
 initCobraToolbox(false);
 changeCobraSolver('gurobi', 'all');
%load('D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/environment')
%restoreEnvironment(environment);

fileName = 'Ecoli_sink_no_media_aerobic_iML1515.mat';
TM_0 = readCbModel(fileName);

% fileName = 'TR_based_modified_ecoli.mat';
% TM_0 = readCbModel(fileName);

TM_0.c(:) = 0;

cd(curr_wd)
TM_0.lb(181) = -readvars("Exch_G.csv");   %iML1515
TM_0.lb(1982) = -readvars("Exch_O.csv");

TM_0.ub(2713:2778) = 10; % iML1515

[min_x_model_round_0, max_x_model_round_0] = fluxVariability(TM_0);

rxn_abbrev = TM_0.rxns;
minimum_flux_round_0 = num2cell(min_x_model_round_0);
maximum_flux_round_0 = num2cell(max_x_model_round_0);

fva_op_round_0 = [rxn_abbrev,minimum_flux_round_0, maximum_flux_round_0];

cd(curr_wd)
writecell(fva_op_round_0,"FVA_1b_obj_0_P1.xlsx")

% part - 2

TM_1 = TM_0;

TM_1.c(2669)=1; % for iML1515


[min_x_model_round_1, max_x_model_round_1] = fluxVariability(TM_1);

rxn_abbrev = TM_1.rxns;
minimum_flux_round_1 = num2cell(min_x_model_round_1);
maximum_flux_round_1 = num2cell(max_x_model_round_1);

fva_op_round_1 = [rxn_abbrev,minimum_flux_round_1, maximum_flux_round_1];

cd(curr_wd)
writecell(fva_op_round_1,"FVA_to_check_P1.xlsx")

TM_1.c(2669)=1; % for iML1515

sol = optimizeCbModel(TM_1);
solv = num2cell(sol.v);
fba_sol = [rxn_abbrev,solv];
writecell(fba_sol,"FBA_to_check_P1.csv")









