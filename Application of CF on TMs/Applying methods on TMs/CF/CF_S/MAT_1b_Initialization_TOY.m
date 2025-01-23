% matlab function - Initialization
% part - 1 
%addpath('/data_birdshire/cobratoolbox');
%addpath('/opt/gurobi1200/linux64/matlab');

curr_wd = 'D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_S_1/';
cd(curr_wd)


% some pre-requisites 
 initCobraToolbox(false);
 changeCobraSolver('gurobi', 'all');

%load('D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/environment')
%restoreEnvironment(environment);

fileName = 'TM_0_obj.mat'; 
TM_0 = readCbModel(fileName);

eb_var = readvars("Exch_b.csv");

TM_0.lb(1)=0;
TM_0.ub(1)=eb_var;

TM_0.lb(2)=0;
TM_0.lb(3)=0;
TM_0.lb(4)=0;
TM_0.lb(5)=0;
TM_0.lb(6)=0;
TM_0.lb(7)=0;
TM_0.lb(8)=0;

TM_0.ub(2)=10000;
TM_0.ub(3)=10000;
TM_0.ub(4)=10000;
TM_0.ub(5)=10000;
TM_0.ub(6)=10000;
TM_0.ub(7)=10000;
TM_0.ub(8)=10000;


[min_x_model_round_0, max_x_model_round_0] = fluxVariability(TM_0);

rxn_abbrev = TM_0.rxns;
minimum_flux_round_0 = num2cell(min_x_model_round_0);
maximum_flux_round_0 = num2cell(max_x_model_round_0);

fva_op_round_0 = [rxn_abbrev,minimum_flux_round_0, maximum_flux_round_0];

cd(curr_wd)
writecell(fva_op_round_0,"FVA_1b_obj_0.xlsx")

% part - 2

TM_1 = TM_0;
TM_1.lb(1)=0;
TM_1.ub(1)=eb_var;

TM_0.lb(2)=0;
TM_0.lb(3)=0;
TM_0.lb(4)=0;
TM_0.lb(5)=0;
TM_0.lb(6)=0;
TM_0.lb(7)=0;
TM_0.lb(8)=0;


TM_0.ub(2)=10000;
TM_0.ub(3)=10000;
TM_0.ub(4)=10000;
TM_0.ub(5)=10000;
TM_0.ub(6)=10000;
TM_0.ub(7)=10000;
TM_0.ub(8)=10000;


 TM_1.c(4)=1;

[min_x_model_round_1, max_x_model_round_1] = fluxVariability(TM_1);

rxn_abbrev = TM_1.rxns;
minimum_flux_round_1 = num2cell(min_x_model_round_1);
maximum_flux_round_1 = num2cell(max_x_model_round_1);

fva_op_round_1 = [rxn_abbrev,minimum_flux_round_1, maximum_flux_round_1];

cd(curr_wd)
writecell(fva_op_round_1,"FVA_to_check.xlsx")

TM_1.c(4)=1;

sol = optimizeCbModel(TM_1);
writematrix(sol.v,"FBA_to_check.csv")

%%%% Final FBA --- only for naive case - WT, C and D KO
%  TM_1.c(4)=1;
%  sol = optimizeCbModel(TM_1);
%  writematrix(sol.v,"FBA_to_check.csv")



