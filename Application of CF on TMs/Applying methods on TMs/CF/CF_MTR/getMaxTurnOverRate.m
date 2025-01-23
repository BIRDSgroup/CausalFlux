function maxTO = getMaxTurnOverRate(model,mets)
% Gets the maximum possible turnover rate for the metabolites in the given model 
%
% USAGE:
%              maxTO = getMaxTurnOverRate(model,mets)
%
% INPUTS:
%              model: COBRA model structure
%              mets:  metIDs for which the maximum turnover rate has to be
%                     found
%
% OUTPUT:
%              maxTO:  maximum turnover rate for all the metablites in
%                      'mets'
maxTO = zeros(numel(mets),1); % zero array for storing the maximum turnover rates
sol = optimizeCbModel(model);
tempModel = model;
Vbio = sol.f;
 tempModel.lb(find(tempModel.c))=Vbio;
% tempModel.lb(find(tempModel.c))=0;
[modelIrrev, ~,~,irrev2rev] = convertToIrreversible(model);
for m=1:numel(mets)
    met = mets{m};
    met_id = find(ismember(modelIrrev.mets,met));
    maxTO(m) = 0.5*maximizeTurnOver(modelIrrev,met_id,irrev2rev);
end
end


function TO = maximizeTurnOver(model,met_id,irrev2rev)
[m,n] = size(model.S);
MetRxnsID = find(model.S(met_id,:)~=0);
nMetRxns = nnz(model.S(met_id,:));
% objective
f = [abs(model.S(met_id,:)');zeros(nMetRxns,1)];

% equalities
Aeq = [model.S, sparse(m,nMetRxns)];
beq = zeros(m,1);
csenseeq = repmat('E',m,1); % equality

% inequalities
temp1 = zeros(n,1);
temp1(MetRxnsID) =1;
temp1 = spdiag(temp1);
temp1(sum(temp1,2)==0,:)=[];
temp2 = -1*spdiag(model.ub(MetRxnsID));
Aineq1 = [temp1,temp2];
bineq1 = zeros(numel(MetRxnsID),1);
csenseineq1 = repmat('L',numel(MetRxnsID),1); % lesser than


Aineq2=[];bineq2=[];csenseineq2=[];
matchRevTemp = irrev2rev(MetRxnsID);
U_matchRevTemp = unique(matchRevTemp);
for i=1:numel(U_matchRevTemp)
    temp_ = U_matchRevTemp(i);
    temp1 = sparse(1,n);
    temp2 = sparse(1,nMetRxns); 
    temp2(matchRevTemp==temp_)=1;
    Aineq2 = [Aineq2;temp1,temp2];
    bineq2 = [bineq2;1];
    csenseineq2 = [csenseineq2;'L']; % lesser than
end

% bounds
lb = [model.lb;zeros(nMetRxns,1)];
ub = [model.ub;ones(nMetRxns,1)];

% Set up MILP problem
MILPproblem.A=[Aeq;Aineq1;Aineq2];
MILPproblem.b=[beq;bineq1;bineq2];
MILPproblem.lb=lb;
MILPproblem.ub=ub;
MILPproblem.c=full(f);
MILPproblem.osense=-1;%maximise
MILPproblem.vartype = [repmat('C',n,1);repmat('B',nMetRxns,1)];
MILPproblem.csense = [csenseeq; csenseineq1; csenseineq2];
solution = solveCobraMILP(MILPproblem);
stat=solution.stat;
if stat~=1
    TO = [];
else
    TO = solution.obj;
end
end