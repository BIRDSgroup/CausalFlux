% matlab function - step 2d  - part 1
% initCobraToolbox(false);
% changeCobraSolver('gurobi', 'all');

load('D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/environment')
restoreEnvironment(environment);

curr_wd = 'D:\work\Integrated_network_model\Ecoli_intg_ntwk\metabolic_aspect\Auto_RUN\Causal_Surgery\Parallel_Runs\CF_S';
cd(curr_wd)
cp_r1 = "CP_round_P1_i.xlsx";
cp_r1 = readtable(cp_r1);

gene_cp_r1.value = cp_r1.Probability ; gene_cp_r1.gene = cp_r1.Bigg_symb;


% fileName = 'Ecoli_model_aerobic_sink_no_media.mat'; 
%  fileName = 'TR_based_modified_ecoli.mat';

fileName = 'Ecoli_sink_no_media_aerobic_iML1515.mat';
TM_0 = readCbModel(fileName);

% % In case of aerobic condition (higher glucose uptake)
%  TM_0.lb(975) = -21.2;
%  TM_0.lb(181) = -10;   %iML1515

% %% In case of Anaerobic condition
% % O2 exchange
%  TM_0.lb(864) = 0;
% % % Glucose 
%  TM_0.lb(975) = -21.2;

% %% When using TR EcoMac MM
% % O2 exchange
% TM_0.lb(933) = 0;
% % % Glucose 
% TM_0.lb(849) = -21.2;

% % % %% In case of Indole flux predictions
% % % O2 exchange
% TM_0.lb(864) = -13;
% % % % Glucose
% TM_0.lb(975) = -9.5;

 cd('D:\work\Integrated_network_model\Ecoli_intg_ntwk\metabolic_aspect\Auto_RUN\Causal_Surgery\CF_S\MM\LB_media\')
 LBmediaconstraints = readtable("LB_media_constraints_iML1515.csv");

 TM_0.ub(LBmediaconstraints.media_list)=-LBmediaconstraints.media_list_ub;
% 
% TM_0.lb(864) = -14.5;

%  TM_0.lb(181) = -8.5;   %iML1515
%  TM_0.lb(1982) = -14.5;

% TM_0.lb(181) = -2.7;   %iML1515
%  TM_0.lb(1982) = -6.71;

 cd(curr_wd)
TM_0.lb(181) = -readvars("Exch_G.csv");   %iML1515
TM_0.lb(1982) = -readvars("Exch_O.csv");

%  TM_0.ub(2383:2448) = 1000;


TM_0.ub(2713:2778) = 10;  % iML1515

gpr_eval_r1 = mapExpressionToReactions(TM_0,gene_cp_r1);
%[expressionRxns ,parsedGPR, gene_used]= mapExpressionToReactions(TM_0,gene_cp_r1);
TF = isnan(gpr_eval_r1);
gpr_eval_r1(TF)=1;

cd(curr_wd)
writematrix(gpr_eval_r1,'GPR_eval_round_P1_i.xlsx')


