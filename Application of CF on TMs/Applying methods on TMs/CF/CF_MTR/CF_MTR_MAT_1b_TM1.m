%% Initialization for MT approach 
%addpath('/data_birdshire/cobratoolbox');
%addpath('/opt/gurobi1200/linux64/matlab');

curr_wd = 'D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_MTR_1/';
cd(curr_wd)

 initCobraToolbox(false);
 changeCobraSolver('gurobi', 'all');

%load('D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/environment')
%restoreEnvironment(environment);

fileName = 'TM_0_obj.mat'; 
mm = readCbModel(fileName);


% to add sink metabolites
extra_mets = {
             'm2[c]',...
             'm4[c]',...
             'm5[c]',...
};


Ml = extra_mets;

e = find(ismember(mm.mets,Ml)); 

mm.c(4)=1;

eb_var = readvars("Exch_b.csv");

mm.lb(1)=0;
mm.ub(1)=eb_var;

mm.lb(2)=0;
mm.lb(3)=0;
mm.lb(4)=0;
mm.lb(5)=0;
mm.lb(6)=0;
mm.lb(7)=0;
mm.lb(8)=0;

mm.ub(2)=10000;
mm.ub(3)=10000;
mm.ub(4)=10000;
mm.ub(5)=10000;
mm.ub(6)=10000;
mm.ub(7)=10000;
mm.ub(8)=10000;


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
writecell(fva_op_round_0,"FVA_1b_obj_0.xlsx")

mm_1 = mm;

mm_1.c(4)=1;

[min_x_model_round_1, max_x_model_round_1] = fluxVariability(mm_1);

rxn_abbrev = mm.rxns;
minimum_flux_round_1 = num2cell(min_x_model_round_1);
maximum_flux_round_1 = num2cell(max_x_model_round_1);

fva_op_round_1 = [rxn_abbrev,minimum_flux_round_1, maximum_flux_round_1];

cd(curr_wd)
writecell(fva_op_round_1,"FVA_to_check.xlsx")

mm_1.c(4)=1;

sol = optimizeCbModel(mm_1);
writematrix(sol.v,"FBA_to_check.csv")



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








