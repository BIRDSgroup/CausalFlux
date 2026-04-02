curr_wd = 'D:\work\Integrated_network_model\Git_hub_codes\Ecoli\CF';
cd(curr_wd)
initCobraToolbox(false);
changeCobraSolver('gurobi', 'all');

cp_r1 = "CP_round_P1_i.xlsx";
cp_r1 = readtable(cp_r1);

gene_cp_r1.value = cp_r1.Probability ; gene_cp_r1.gene = cp_r1.Bigg_symb;

fileName = 'NEW_iML1515_modified_sink.mat';
TM_0 = readCbModel(fileName);

cd(curr_wd)
exch_rxns_dt = readtable("Exch_iml1515.csv");

TM_0.lb(exch_rxns_dt.Var1) = 0;

cd(curr_wd)
essen_exch_rxns_dt = readtable("essential_exch_rxns_ids_iml1515.csv");
TM_0.lb(essen_exch_rxns_dt.Var1)= -0.001;


cd(curr_wd)
LBmediaconstraints = readtable("LB_media_constraints_iML1515.csv");
% LBmediaconstraints = readtable("M9_serine_galacturonate.csv");
% LBmediaconstraints = readtable("TSBmed.csv");

 TM_0.ub(LBmediaconstraints.media_list)=-LBmediaconstraints.media_list_ub;


cd(curr_wd)
TM_0.lb(181) = -readvars("Exch_G.csv");   %iML1515 glucose exchange reaction 
TM_0.lb(1982) = -readvars("Exch_O.csv");  %iML1515 oxygen exchange reaction 

TM_0.ub(2713:2731) = 10; % iML1515 extra sink reactions added based on metabolic feedback information 

gpr_eval_r1 = mapExpressionToReactions(TM_0,gene_cp_r1);
TF = isnan(gpr_eval_r1);
gpr_eval_r1(TF)=1;

cd(curr_wd)
writematrix(gpr_eval_r1,'GPR_eval_round_P1_i.xlsx')


