%% Initialization for MT approach 

curr_wd = 'D:\work\Integrated_network_model\Ecoli_intg_ntwk\metabolic_aspect\Auto_RUN\Causal_Surgery\Parallel_Runs\CF_MTR';
cd(curr_wd)

 initCobraToolbox(false);
 changeCobraSolver('gurobi', 'all');
%load('D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/environment')
%restoreEnvironment(environment);

% fileName = 'Ecoli_aerobic_no_sink.mat';  
fileName = 'Ecoli_no_sink_no_media_aerobic_iML1515.mat';

mm = readCbModel(fileName);

% % In case of aerobic condition (higher glucose uptake)
%  mm.lb(975) = -21.2;
% mm.lb(975) = -14;

% % In case of Anaerobic condition
% % O2 exchange
% mm.lb(864) = 0;
% % Glucose 
% mm.lb(975) = -21.2;

% to add sink metabolites
extra_mets = {'cytd_c',...
             'fe2_c',...
             'arg__L_c',...
             'leu__L_c',...
             'lys__L_c',...
             'thym_c',...
             'ura_c',...
             'camp_c',...
             'fdp_c',...
             'f1p_c',...
             'ade_c',...
             'for_c',...
             'glcn_c',...
             'gly_c',...
             'glx_c',...
             'gua_c',...
             'acgam6p_c',...
             'acser_c',...
             'ametam_c',...
             'alltn_c',...
             'h2s_c',...
             'tsul_c',...
             '2mcit_c',...
             'alac__S_c',...
             'lac__L_c',...
             '2ahbut_c',...
             '2dr5p_c',...
             'dhpppn_c',...
             'pppn_c',...
             'dhptd_c',...
             '5dglcn_c',...
             'ac_c',...
             'aps_c',...
             'chol_c',...
             'cynt_c',...
             'fruur_c',...
             'rib__D_c',...
             'ser__D_c',...
             'xyl__D_c',...
             'glyc_c',....
             'glyclt_c',...
             'hxan_c',...
             'idon__L_c',...
             'ala__L_c',...
             'arab__L_c',...
             'asn__L_c',...
             'crn__D_c',...
             'fc1p_c',...
             'hcys__L_c',...
             'lcts_c',...
             'phe__L_c',...
             'rmn_c',...
             'trp__L_c',...
             'tyr__L_c',...
             'malttr_c',...
             'melib_c',...
             'acnam_c',...
             'na1_c',...
             'phaccoa_c',...
             'pep_c',...
             'pyr_c',...
             'tre_c',...
             'tre6p_c',...
             'xtsn_c',...
             'gal_c',...
             'glyc3p_c',...
};


Ml = extra_mets;

e = find(ismember(mm.mets,Ml)); 

% tr = maxTResti(mm,Ml);

% mm.c(926)=1;



 cd('D:\work\Integrated_network_model\Ecoli_intg_ntwk\metabolic_aspect\Auto_RUN\Causal_Surgery\Ishii_dataset\MM\LB_media\')
 LBmediaconstraints = readtable("LB_media_constraints_iML1515.csv");

mm.lb(LBmediaconstraints.media_list)=-LBmediaconstraints.media_list_ub;

cd(curr_wd)

mm.c(2669)=1;


cd(curr_wd)
mm.lb(181) = -readvars("Exch_G.csv");   %iML1515
mm.lb(1982) = -readvars("Exch_O.csv");

% mm.lb(181) = -3.28;   %iML1515
%  mm.lb(1982) = -1.06;

tr = getMaxTurnOverRate(mm,Ml);




% tr_ = transpose(tr);
en = mm.metNames(e);
es = mm.mets(e);
% tr_1 = num2cell(tr_);
% C2 = [num2cell(e),es,tr_1];
tr_1 = num2cell(tr);
C2 = [num2cell(e),es,tr_1];

cd(curr_wd)
writecell(C2,"Initial_Max_Turnover.csv")

% initCobraToolbox(false);
% changeCobraSolver('gurobi', 'all');


mm.c(:) = 0;

[min_x_model_round_0, max_x_model_round_0] = fluxVariability(mm);

rxn_abbrev = mm.rxns;
minimum_flux_round_0 = num2cell(min_x_model_round_0);
maximum_flux_round_0 = num2cell(max_x_model_round_0);

fva_op_round_0 = [rxn_abbrev,minimum_flux_round_0, maximum_flux_round_0];

cd(curr_wd)
writecell(fva_op_round_0,"FVA_1b_obj_0_P2.xlsx")

mm_1 = mm;

% mm_1.c(926)=1;  % for our IAF1260 - Biomass
%  TM_1.c(1005)=1;  % for TR Ecoli - Biomass

% TM_1.c(2267)=1;  % for our IAF1260 - Indole
% TM_1.c(2272)=1;  % for TR Ecoli - Biomass

mm_1.c(2669)=1;

% [min_x_model_round_1, max_x_model_round_1] = fluxVariability(mm_1);
% 
% rxn_abbrev = mm.rxns;
% minimum_flux_round_1 = num2cell(min_x_model_round_1);
% maximum_flux_round_1 = num2cell(max_x_model_round_1);
% 
% fva_op_round_1 = [rxn_abbrev,minimum_flux_round_1, maximum_flux_round_1];
% 
% cd(curr_wd)
% writecell(fva_op_round_1,"FVA_to_check_P2.xlsx")




% mm_1.c(926)=1;  % for our IAF1260 - Biomass
%  TM_1.c(1005)=1;  % for TR Ecoli - Biomass

% TM_1.c(2267)=1;  % for our IAF1260 - Indole
% TM_1.c(2272)=1;  % for TR Ecoli - Biomass

mm_1.c(2669)=1;

%TM_1.ub(2383:2448) = 1000;


sol = optimizeCbModel(mm_1);
writematrix(sol.v,"FBA_to_check_P2.csv")



mvec = find(ismember(mm_1.mets,Ml)); 
N = length(mvec);
iter_TR = zeros(1,N); 

for i = 1:length(mvec)
    
   iter_TR(i) = 0.5 * (sum(abs(mm_1.S(mvec(i),:).*transpose(sol.v))));   %  resolve this issue
    
      
end

iter_TR_ = transpose(iter_TR);

% iter_TR_ = getMaxTurnOverRate(mm_1,Ml);
en = mm.metNames(mvec);
es = mm.mets(mvec);
iter_TR_1 = num2cell(iter_TR_);
C4 = [num2cell(mvec),es,iter_TR_1];

cd(curr_wd)
writecell(C4,"Iteration_Turnover.csv")








