%% Iteration step - 1 for MT approach 

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
% mm.lb(975) = -21.2;
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


% mm.c(926) = 1;
% sol = optimizeCbModel(mm);



 cd('D:\work\Integrated_network_model\Ecoli_intg_ntwk\metabolic_aspect\Auto_RUN\Causal_Surgery\Ishii_dataset\MM\LB_media\')
 LBmediaconstraints = readtable("LB_media_constraints_iML1515.csv");

  mm.lb(LBmediaconstraints.media_list)=-LBmediaconstraints.media_list_ub;

cd(curr_wd)
 mm.lb(181) = -readvars("Exch_G.csv");   %iML1515
 mm.lb(1982) = -readvars("Exch_O.csv");

% mm.lb(181) = -2.86;   %iML1515
%  mm.lb(1982) = -6.64;

curr_wd = 'D:\work\Integrated_network_model\Ecoli_intg_ntwk\metabolic_aspect\Auto_RUN\Causal_Surgery\Parallel_Runs\CF_MTR';
cd(curr_wd)

Updated_FVA_round_i = "Updated_FVA_round_P2_i.xlsx";
Updated_FVA_round_i = readtable(Updated_FVA_round_i, "VariableNamingRule","preserve");

x_model_round_i = mm;

% x_model_round_i.c(926)=1; % IAF1260 - BIOMASS
% x_model_round_i.c(1005)=1; % TR EcoMac - Biomass

% x_model_round_i.c(2267)=1;  % for our IAF1260 - Indole
% x_model_round_i.c(2272)=1;  % for TR Ecoli - Indole

x_model_round_i.c(2669)=1;

x_model_round_i.ub = Updated_FVA_round_i.new_upper_bounds;
x_model_round_i.lb = Updated_FVA_round_i.new_lower_bounds;

% [min_x_model_ri, max_x_model_ri] = fluxVariability(x_model_round_i); %with
% %round 0 FVA used for upper bounds evaluation
% 
% %[min_x_model_r2_, max_x_model_r2_] = fluxVariability(x_model_round_2);
% 
% rxn_abbrev = x_model_round_i.rxns;
% minimum_flux_ri = num2cell(min_x_model_ri);
% maximum_flux_ri = num2cell(max_x_model_ri);
% fva_op_round_i = [rxn_abbrev,minimum_flux_ri, maximum_flux_ri]; %with
% 
% 
% cd(curr_wd)
% writecell(fva_op_round_i,"FVA_to_check_P2.xlsx")



% x_model_round_i.c(926)= 1;  % IAF1260 - BIOMASS
% x_model_round_i.c(1005)=1; % TR EcoMac - Biomass
% x_model_round_i.c(2267)=1;  % for our IAF1260 - Indole
% x_model_round_i.c(2272)=1;  % for TR Ecoli - Indole

x_model_round_i.c(2669)= 1;

sol = optimizeCbModel(x_model_round_i);

if sol.stat==1
 writematrix(sol.v,"FBA_to_check_P2.csv")
 
 
 mvec = find(ismember(x_model_round_i.mets,Ml)); 
 N = length(mvec);
 Curr_TR = zeros(1,N); 
 
 % step 2 - For loop for each mets present in the Mets_L
 for i = 1:length(mvec)
     
    Curr_TR(i) = 0.5 * (sum(abs(mm.S(mvec(i),:).*transpose(sol.v))));   %  resolve this issue
     
       
 end
 
 
 Curr_TR_ = transpose(Curr_TR);
 % Curr_TR_ = getMaxTurnOverRate(x_model_round_i,Ml);
 
 e = find(ismember(x_model_round_i.mets,Ml));
 es = x_model_round_i.mets(mvec);
 tr_1 = num2cell(Curr_TR_);
 C3 = [num2cell(e),es,tr_1];
 
 cd(curr_wd)
 writecell(C3,"Iteration_Turnover.csv")
end



















