% matlab function - Initialization
% part - 1 
curr_wd = 'D:\work\Integrated_network_model\Toy_model\auto_res';
cd(curr_wd)


% some pre-requisites 
initCobraToolbox(false);
changeCobraSolver('gurobi', 'all');

fileName = 'TM_0_obj.mat'; 
TM_0 = readCbModel(fileName);


[min_x_model_round_0, max_x_model_round_0] = fluxVariability(TM_0);

rxn_abbrev = TM_0.rxns;
minimum_flux_round_0 = num2cell(min_x_model_round_0);
maximum_flux_round_0 = num2cell(max_x_model_round_0);

fva_op_round_0 = [rxn_abbrev,minimum_flux_round_0, maximum_flux_round_0];

cd(curr_wd)
writecell(fva_op_round_0,"FVA_1b_obj_0.xlsx")

% part - 2

TM_1 = TM_0;
TM_1.c(4)=1;

[min_x_model_round_1, max_x_model_round_1] = fluxVariability(TM_1);

rxn_abbrev = TM_1.rxns;
minimum_flux_round_1 = num2cell(min_x_model_round_1);
maximum_flux_round_1 = num2cell(max_x_model_round_1);

fva_op_round_1 = [rxn_abbrev,minimum_flux_round_1, maximum_flux_round_1];

cd(curr_wd)
writecell(fva_op_round_1,"FVA_to_check.xlsx")
