% matlab function - step 2d  - part 1
%addpath('/data_birdshire/cobratoolbox');
%addpath('/opt/gurobi1200/linux64/matlab');


 initCobraToolbox(false);
 changeCobraSolver('gurobi', 'all');

%load('D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/environment')
%restoreEnvironment(environment);

curr_wd = 'D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_S_1/';
cd(curr_wd)
cp_r1 = "CP_round_i.xlsx";
cp_r1 = readtable(cp_r1);

gene_cp_r1.value = cp_r1.Probability ; gene_cp_r1.gene = cp_r1.Bigg_symb;


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

gpr_eval_r1 = mapExpressionToReactions(TM_0,gene_cp_r1);
%[expressionRxns ,parsedGPR, gene_used]= mapExpressionToReactions(TM_0,gene_cp_r1);
TF = isnan(gpr_eval_r1);
gpr_eval_r1(TF)=1;

cd(curr_wd)
writematrix(gpr_eval_r1,'GPR_eval_round_i.xlsx')
