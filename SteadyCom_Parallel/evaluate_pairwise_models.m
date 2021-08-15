close all ; clc; clear
tStart = tic ; 
models_55 = load('models55_SteadyCom.mat') ;
mediumMets = models_55.mediumMets ;

final_pairwise_results = cell(1,1) ; 
fields_name = cell(1,1) ;

for l = 2:2
    
    comm_model = models_55.comm_models{l-1} ;
    nameList = comm_model.modelID ;
    model_1 = models_55.all_models{1} ; 
    model_2 = models_55.all_models{l} ; 
    models = {model_1 ; model_2} ; 
    name_model = strcat('Comm_model_',num2str(l),'__',comm_model.description) ; 
    final_results = run_KO_test_pairs(models,nameList, mediumMets);
    final_pairwise_results{1} = final_results ; 
    fields_name{1} = name_model ;
    name_struct = strcat(name_model,'.mat') ; 
    save(name_struct,'final_pairwise_results')
end
%evaluate_pairwise_result = cell2struct(final_pairwise_results, fields_name, 1) 
tEnd = toc(tStart) 
    