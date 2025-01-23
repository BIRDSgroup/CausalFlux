% matlab function - step 2a  - FVA with 0 objective function 
%addpath('/data_birdshire/cobratoolbox');
%addpath('/opt/gurobi1200/linux64/matlab');

curr_wd = 'D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_S_1/';
cd(curr_wd)

Updated_FVA_round_i = "Updated_FVA_round_i.xlsx";
Updated_FVA_round_i = readtable(Updated_FVA_round_i, "VariableNamingRule","preserve");

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


x_model_round_i = TM_0;
x_model_round_i.c(4)=1;
x_model_round_i.ub = Updated_FVA_round_i.new_upper_bounds;
[min_x_model_ri, max_x_model_ri] = fluxVariability(x_model_round_i); %with
%round 0 FVA used for upper bounds evaluation

%[min_x_model_r2_, max_x_model_r2_] = fluxVariability(x_model_round_2);

rxn_abbrev = x_model_round_i.rxns;
minimum_flux_ri = num2cell(min_x_model_ri);
maximum_flux_ri = num2cell(max_x_model_ri);
fva_op_round_i = [rxn_abbrev,minimum_flux_ri, maximum_flux_ri]; %with


cd(curr_wd)
writecell(fva_op_round_i,"FVA_to_check.xlsx")

x_model_round_i.c(4)=1;
sol = optimizeCbModel(x_model_round_i);
writematrix(sol.v,"FBA_to_check.csv")

x_model_round_i.c(4)=1;

sol = optimizeCbModel(x_model_round_i);
writematrix(sol.v,"FBA_to_check.csv")


