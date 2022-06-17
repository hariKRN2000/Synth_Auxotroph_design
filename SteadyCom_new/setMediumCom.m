function newModel = setMediumCom(model,mediumMets)

model_old = model;
mediumRxns = intersect(model.indCom.EXcom,find(ismember(model.rxns,findRxnsFromMets(model,mediumMets))));

%Set UB for all mets
model.ub(model.indCom.EXcom) = 1000; % Unconstrain upper bounds on all exchange reactions except noSecreteMets
noSecreteMets = {'fe2_e[u]','fe3_e[u]','fe_e[u]','n2_e[u]'};
model.ub(intersect(model.indCom.EXcom,find(ismember(model.rxns,findRxnsFromMets(model,noSecreteMets))))) = 0;

%Set LB for all mets
model.lb(model.indCom.EXcom) = 0;
for i = 1:length(mediumRxns)
    model.lb(mediumRxns(i)) = -2;
end

model.lb(intersect(model.indCom.EXcom,find(ismember(model.rxns,findRxnsFromMets(model,{'glc__D_e[u]'}))))) = -10;
model.lb(intersect(model.indCom.EXcom,find(ismember(model.rxns,findRxnsFromMets(model,{'malt_e[u]'}))))) = -0;
newModel = model;
end
