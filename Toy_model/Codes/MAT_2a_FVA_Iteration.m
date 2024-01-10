% matlab function - step 2a  - FVA with 0 objective function 
curr_wd = 'D:\work\Integrated_network_model\Toy_model\auto_res';
cd(curr_wd)

Updated_FVA_round_i = "Updated_FVA_round_i.xlsx";
Updated_FVA_round_i = readtable(Updated_FVA_round_i, "VariableNamingRule","preserve");

% some pre-requisites 
initCobraToolbox(false);
changeCobraSolver('gurobi', 'all');

fileName = 'TM_0_obj.mat'; 
TM_0 = readCbModel(fileName);




x_model_round_i = TM_0;
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
