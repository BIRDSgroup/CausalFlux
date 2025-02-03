% matlab function - step 2a  - FVA with 0 objective function 
curr_wd = 'D:\work\Integrated_network_model\Ecoli_intg_ntwk\metabolic_aspect\Auto_RUN\Causal_Surgery\Parallel_Runs\CF_S';
cd(curr_wd)

Updated_FVA_round_i = "Updated_FVA_round_P1_i.xlsx";
Updated_FVA_round_i = readtable(Updated_FVA_round_i, "VariableNamingRule","preserve");

% some pre-requisites 
 initCobraToolbox(false);
 changeCobraSolver('gurobi', 'all');
%load('D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/environment')
%restoreEnvironment(environment);


% fileName = 'Ecoli_model_aerobic_sink_no_media.mat'; 
% fileName = 'TR_based_modified_ecoli.mat';

fileName = 'Ecoli_sink_no_media_aerobic_iML1515.mat';

TM_0 = readCbModel(fileName);

 cd('D:\work\Integrated_network_model\Ecoli_intg_ntwk\metabolic_aspect\Auto_RUN\Causal_Surgery\CF_S\MM\LB_media\')
 LBmediaconstraints = readtable("LB_media_constraints_iML1515.csv");

 TM_0.ub(LBmediaconstraints.media_list)=-LBmediaconstraints.media_list_ub;

 cd(curr_wd)
TM_0.lb(181) = -readvars("Exch_G.csv");   %iML1515
TM_0.lb(1982) = -readvars("Exch_O.csv");

TM_0.ub(2713:2778) = 10;  % iML1515

x_model_round_i = TM_0;

x_model_round_i.c(2669)=1; % for iML1515

x_model_round_i.ub = Updated_FVA_round_i.new_upper_bounds;
x_model_round_i.lb = Updated_FVA_round_i.new_lower_bounds;

[min_x_model_ri, max_x_model_ri] = fluxVariability(x_model_round_i); %with
%round 0 FVA used for upper bounds evaluation

%[min_x_model_r2_, max_x_model_r2_] = fluxVariability(x_model_round_2);

rxn_abbrev = x_model_round_i.rxns;
minimum_flux_ri = num2cell(min_x_model_ri);
maximum_flux_ri = num2cell(max_x_model_ri);
fva_op_round_i = [rxn_abbrev,minimum_flux_ri, maximum_flux_ri]; %with


cd(curr_wd)
writecell(fva_op_round_i,"FVA_to_check_P1.xlsx")
 
x_model_round_i.c(2669)=1; % for iML1515

sol = optimizeCbModel(x_model_round_i);
solv = num2cell(sol.v);
fba_sol = [rxn_abbrev, solv];
writecell(fba_sol,"FBA_to_check_P1.csv")




