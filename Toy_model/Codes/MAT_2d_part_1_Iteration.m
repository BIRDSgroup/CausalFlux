% matlab function - step 2d  - part 1
initCobraToolbox(false);
changeCobraSolver('gurobi', 'all');

curr_wd = 'D:\work\Integrated_network_model\Toy_model\auto_res';
cd(curr_wd)
cp_r1 = "CP_round_i.xlsx";
cp_r1 = readtable(cp_r1);

gene_cp_r1.value = cp_r1.Probability ; gene_cp_r1.gene = cp_r1.Bigg_symb;


fileName = 'TM_0_obj.mat'; 
TM_0 = readCbModel(fileName);

gpr_eval_r1 = mapExpressionToReactions(TM_0,gene_cp_r1);
%[expressionRxns ,parsedGPR, gene_used]= mapExpressionToReactions(TM_0,gene_cp_r1);
TF = isnan(gpr_eval_r1);
gpr_eval_r1(TF)=1;

cd(curr_wd)
writematrix(gpr_eval_r1,'GPR_eval_round_i.xlsx')
