clear; clc; 
ecoli = readCbModel('iML1515.mat');
yeast = readCbModel('iML1515.mat');

%minimal medium
mediumMets = {'his__L_e[u]';'ala__L_e[u]';'leu__L_e[u]';'met__L_e[u]';'gln__L_e[u]';'ile__L_e[u]';'arg__L_e[u]';'trp__L_e[u]';'val__L_e[u]';'4abz_e[u]'; 'btn_e[u]'; 'ca2_e[u]';'cbl1_e[u]';'cl_e[u]';'cobalt2_e[u]';'cu2_e[u]';'fe2_e[u]';'fe3_e[u]';'fol_e[u]';'h2_e[u]';'h2o_e[u]';'k_e[u]';'mg2_e[u]';'mn2_e[u]';'mobd_e[u]';'na1_e[u]';'ncam_e[u]';'nh4_e[u]';'ni2_e[u]';'no3_e[u]';'o2_e[u]';'pi_e[u]';'pnto__R_e[u]';'ribflv_e[u]';'so4_e[u]';'zn2_e[u]';'nac_e[u]';'sel_e[u]';'adocbl_e[u]';'pb2_e[u]';'co2_e[u]';'so3_e[u]';'h2s_e[u]';'meoh_e[u]';'pime_e[u]';'pheme_e[u]';'gua_e[u]';'orot_e[u]';'xan_e[u]';'thymd_e[u]';'ura_e[u]';'thm_e[u]'};
 
all_models = {ecoli, yeast };
  all_biomass = {char(ecoli.rxns(find(ecoli.c))) ; char(yeast.rxns(find(yeast.c)))} ; % Char array of names of biomass reactions
  modelList = {'ecoli' ,'yeast'};

total_models = nchoosek(length(all_models),2);

comm_models ={};
pairedModelInfo ={};
count=0;
temp = 0;

for i = 1:length(all_models)
    for j = i+1 :length(all_models)
        model1 = all_models{i}; model2 = all_models{j};
        
        exc1 = model1.rxns(findExcRxns(model1)==1);
        exc2 = model2.rxns(findExcRxns(model2)==1);
        model1 = changeRxnBounds(model1,exc1,-1000,'l');
        model2 = changeRxnBounds(model2,exc2,-1000,'l');
        model1 = changeRxnBounds(model1,exc1,1000,'u');
        model2 = changeRxnBounds(model2,exc2,1000,'u');
        models = {model1,model2};
        options.spBm = {all_biomass{i}, all_biomass{j}};        
        options.spAbbr = {modelList{i},modelList{j}};
        options.metExId = ('_e');
        options.sepUtEx = false;
        
        count = count+1;
        ModelCom = createCommModel1(models,options);
        
%for glucose as carbon source
        ModelCom = setMediumCom(ModelCom, mediumMets);
                 
        comm_models{count} = ModelCom;
        % information on joined models
            pairedModelInfo{count, 1} = strcat('pairedModel', '_', modelList{i}, '_', modelList{j}, '.mat');
            pairedModelInfo{count, 2} = modelList{i};
            pairedModelInfo{count, 3} = all_biomass{i};
            pairedModelInfo{count, 4} = modelList{j};
            pairedModelInfo{count, 5} = all_biomass{j};     
    end
end

options = struct();
options.GRguess = 0.1;  % initial guess for max. growth rate
options.GRtol = 1e-6;  % tolerance for final growth rate
options.algorithm = 1;  % use the default algorithm (simple guessing for bounds, followed by matlab fzero)
[sol_WT, result_WT] = SteadyCom(ModelCom);
result_WT.GRmax
