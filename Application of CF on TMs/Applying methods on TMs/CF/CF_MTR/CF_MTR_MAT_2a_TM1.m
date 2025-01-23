%% Iteration step - 1 for MT approach 
%addpath('/data_birdshire/cobratoolbox');
%addpath('/opt/gurobi1200/linux64/matlab');

curr_wd = 'D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/CF_MTR_1/';
cd(curr_wd)

% initCobraToolbox(false);
% changeCobraSolver('gurobi', 'all');

load('D:/work/Integrated_network_model/Toy_model/auto_new_model_current_approach_27_03_24/Causal_Surgery/environment')
restoreEnvironment(environment);

fileName = 'TM_0_obj.mat'; 
mm = readCbModel(fileName);

% to add sink metabolites
extra_mets = {
             'm2[c]',...
             'm4[c]',...
             'm5[c]',...
};

Ml = extra_mets;

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


Updated_FVA_round_i = "Updated_FVA_round_i.xlsx";
Updated_FVA_round_i = readtable(Updated_FVA_round_i, "VariableNamingRule","preserve");

x_model_round_i = mm;

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


x_model_round_i.c(4)= 1;

sol = optimizeCbModel(x_model_round_i);
writematrix(sol.v,"FBA_to_check.csv")


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



















