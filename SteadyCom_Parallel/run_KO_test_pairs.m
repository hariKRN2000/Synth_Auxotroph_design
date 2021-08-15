function final_results = run_KO_test_pairs(models,nameList, mediumMets,commModel)
% If only first inout is given (as a cell), then the function will create a
% CommModel of the given models and evaluate the KOs and give out the
% results. nameList is preffered so that you can name the pairs efficiently
% If third argument is given (list of all medium mets), then the function will evaluate by settin
% the required medium according to the medium mets list
% Note that if you have 4th argument in, 3rd argument is necessary
%close all ;clear; clc; 
tic
% Setting up solver
changeCobraSolver('gurobi')

% Creating Structure to store results
excel_name = 'results_excel.xls' ; 
final_results = struct() ; 
if nargin >= 2  
    final_results(1).Organims = strcat(nameList{1},'_',nameList{2}) ;
    names = strcat([nameList{1}],'_',[nameList{2}])  ;
    excel_name = strcat('results__',names,'__excel.xls') ;
end

% Loading and creating Comm Models
model_1_WT = models{1} ; %readCbModel('iML1515.mat');
model_2_WT = models{2} ; %readCbModel('iMM904.mat') ;
all_models = {model_1_WT,model_2_WT} ; % Creating model list for WT

switch nargin
    case 2
        [comm_models,pairedModelInfo] = create_community(all_models) ; % Creating community model for WT in a set minimal medium
        [sol_WT, result_WT] = SteadyCom(comm_models{1}) ; % Applying SteadyCom to WT community model
        comm_model = comm_models{1} ; 
    case 3
        [comm_models,pairedModelInfo] = create_community(all_models,mediumMets) ; % Creating community model for WT in a defined medium
        [sol_WT, result_WT] = SteadyCom(comm_models{1}) ; % Applying SteadyCom to WT community model
        comm_model = comm_models{1} ;
    case 4
        comm_models = setMediumCom(commModel, mediumMets);  % If input is already a comm model, then just change the medium
        [sol_WT, result_WT] = SteadyCom(comm_models) ; % Applying SteadyCom to WT community model
        comm_model = comm_models ;
end

% Setting medium constraints
%comm_models{1} = changeRxnBounds(comm_models{1},'EX_o2_e(u)',-0,'l'); % Set maximum oxygen uptake
%comm_models{1} = changeRxnBounds(comm_models{1},'EX_glc__D_e(u)',0,'l'); %To remove glucose from medium
%comm_models{1} = changeRxnBounds(comm_models{1},'EX_malt_e(u)',-10,'l'); %To change carbon source

 
GR_WT = result_WT.GRmax ; % Growth rate for WT 
GR_WT_org_1 = result_WT.vBM(1)  % Growth Rate for WT org_1
GR_WT_org_2 = result_WT.vBM(2)  % Growth Rate for WT org_2
wild_type_abundances = result_WT.BM % Abundances

% Storing WT comm model data : 
WT_comm_model = struct() ; 
WT_comm_model.Growth_Rates = result_WT.vBM ; 
WT_comm_model.Abundances = result_WT.BM ; 
final_results(1).WT_Information = WT_comm_model  ;


% Finding Lethal genes for WT org_2
[grRatio_single_org_2, grRateKO_single_org_2] = singleGeneDeletion(model_2_WT);                               
single_lethal_genes_ind_org_2 = find(grRatio_single_org_2 < 0.05 ) ;
% Finding Lethal genes for WT org_1
[grRatio_single_org_1, grRateKO_single_org_1] = singleGeneDeletion(model_1_WT);                               
single_lethal_genes_ind_org_1 = find(grRatio_single_org_1 < 0.05 ) ;
final_results(1).Single_Synthetic_Lethal_deletions_gene_indices.Org_1 = single_lethal_genes_ind_org_1 ;
final_results(1).Single_Synthetic_Lethal_deletions_gene_indices.Org_2 = single_lethal_genes_ind_org_2 ;

% Creating arrays to store important results
result_array_org_1 = zeros(length(single_lethal_genes_ind_org_1),1) ; % array noting which KO grows with WT 
result_array_org_2 = zeros(length(single_lethal_genes_ind_org_2),1) ; % array noting which KO grows with WT 


result_array_Gr_KO_org_1 = zeros(length(single_lethal_genes_ind_org_1),1) ; % ratio of growth rate of KO to WT 
result_array_Gr_KO_org_2 = zeros(length(single_lethal_genes_ind_org_2),1) ; % ratio of growth rate of KO to WT 

result_array_Gr_WT_org_1 = zeros(length(single_lethal_genes_ind_org_2),1) ; % ratio of growth rate of WT with KO to WT 
result_array_Gr_WT_org_2 = zeros(length(single_lethal_genes_ind_org_1),1) ; % ratio of growth rate of WT with KO to WT 

result_array_Ex_org_1 = zeros(length(single_lethal_genes_ind_org_1),1) ; % Export fluxes where KO is org_1
result_array_Ex_org_2 = zeros(length(single_lethal_genes_ind_org_2),1) ; % Export fluxes where KO is org_2

result_array_abundance_org_1 = zeros(length(single_lethal_genes_ind_org_1),1) ; % Abundance of org_1 where KO is org_1
result_array_abundance_org_2 = zeros(length(single_lethal_genes_ind_org_2),1) ; % Abundance of org_1 where KO is org_2

result_array_GRmax_KO_org_1 = zeros(length(single_lethal_genes_ind_org_1),1) ;
result_array_GRmax_KO_org_2 = zeros(length(single_lethal_genes_ind_org_2),1) ;
fault_gene_ind_org_1 = [] ;
fault_gene_ind_org_2 = [] ;
 

environment = getEnvironment();

% Simulating models for KO org_1 : 
parfor i = 1:length(single_lethal_genes_ind_org_1) 
    restoreEnvironment(environment);    
    
    disp(single_lethal_genes_ind_org_1(i))
    disp(i)
    
%      if ismember(i,[41,81]) % &&  ismember(j,error_ind_gene_org_2)
%              %ind = ind + 1 ;
%               continue
%      end
    

    targ_gene_ind = single_lethal_genes_ind_org_1(i) ; 
    ko_gene_org_1 = strcat(model_1_WT.genes(targ_gene_ind) , {'_model_1'}) ;  
    ko_comm_model = deleteModelGenes(comm_model, ko_gene_org_1) ;   
    
    try
    [sol_KO, result_KO] = SteadyCom(ko_comm_model) ; % Applying Steady Come to WT community model
    
    if result_KO.vBM(1) > 0 
        result_array_org_1(i) = 1 ;         
    end
    GR = result_KO.vBM ; Ex = sum(result_KO.Ex) ;
    result_array_Gr_KO_org_1(i) = GR(1)/GR_WT_org_1 ;
    result_array_Gr_WT_org_2(i) = GR(2)/GR_WT_org_2 ;
    result_array_Ex_org_1(i) = Ex ;
    result_array_abundance_org_1(i) = result_KO.BM(1) ; 
    result_array_GRmax_KO_org_1(i) = result_KO.GRmax ;
    
    catch
        fault_gene_ind_org_1 = [fault_gene_ind_org_1; targ_gene_ind] ;
    end
    
end

gene_ind = single_lethal_genes_ind_org_1(find(result_array_org_1)) ; 
gr_KO_org_1 = result_array_Gr_KO_org_1(find(result_array_org_1)) ; 
gr_WT_org_2 = result_array_Gr_WT_org_2(find(result_array_org_1)) ; 
ex_flux = result_array_Ex_org_1(find(result_array_org_1)) ;
grmax = result_array_GRmax_KO_org_1(find(result_array_org_1)) ;
genes = model_1_WT.genes(gene_ind) ; 
abundance_KO = result_array_abundance_org_1(find(result_array_org_1)) ;
Table_org_1 = table(gene_ind,genes, (gr_KO_org_1'*GR_WT_org_1)', (gr_WT_org_2'*GR_WT_org_2)',gr_KO_org_1,grmax,abundance_KO,ex_flux, 'VariableNames',...
             {'Gene Index','Gene', 'KO Growth rate' , 'WT Growth rate', 'Ratio (gr_KO/gr_WT)','GRmax','KO Abundance','Sum of Ex Flux'}) ;
writetable(Table_org_1, excel_name,'Sheet',1) ; 
% Storing information into structures
result_struct_org_1 = table2struct(Table_org_1) ; 
final_results(1).Org_1_KO_Org2_WT_info_table = result_struct_org_1 ;
final_results(1).Org_1_KO_Org2_WT_info.Growth_rate_ratio_Org_1 = result_array_Gr_KO_org_1 ; 
final_results(1).Org_1_KO_Org2_WT_info.Growth_rate_ratio_Org_2 = result_array_Gr_WT_org_2 ; 
final_results(1).Org_1_KO_Org2_WT_info.Exchange_flux_sum = result_array_Ex_org_1 ; 
final_results(1).Org_1_KO_Org2_WT_info.Abundance_Org_1 = result_array_abundance_org_1 ; 

figure
bar(single_lethal_genes_ind_org_1, result_array_Gr_KO_org_1)
xlabel('Index of Gene in org_1 deleted') ; 
ylabel('Growth Rate ratio for Mutant (org_1) in Comm Model')
saveas(gcf, "glucose_org_1_KO_gr_ratio.png");

figure
bar(single_lethal_genes_ind_org_1, result_array_Gr_WT_org_2)
xlabel('Index of Gene in org_1 deleted') ; 
ylabel('Growth Rate ratio for WT (org_2) in Comm Model')
saveas(gcf, "glucose_model_2_WT_gr_ratio.png");


figure
bar(single_lethal_genes_ind_org_1, result_array_Ex_org_1)
xlabel('Index of Gene in org_1 deleted') ; 
ylabel('Exchange Fluxes (mean)')
saveas(gcf, "glucose_org_1_KO_Ex.png");

figure
bar(single_lethal_genes_ind_org_1, result_array_abundance_org_1)
xlabel('Index of Gene in org_1 deleted') ; 
ylabel('org_1 Abundance')
saveas(gcf, "glucose_org_1_KO_Abundance.png");


% Simulating models for KO org_2 : 
parfor i = 1:length(single_lethal_genes_ind_org_2)
    restoreEnvironment(environment);   
    
    targ_gene_ind = single_lethal_genes_ind_org_2(i) ;
    ko_gene_org_2 = strcat(model_2_WT.genes(targ_gene_ind) , {'_model_2'}) ;   
    ko_comm_model = deleteModelGenes(comm_model, ko_gene_org_2) ; 
   
    try
    [sol_KO, result_KO] = SteadyCom(ko_comm_model) ; % Applying Steady Come to WT community model

    
    if result_KO.vBM(1) > 0
        result_array_org_2(i) = 1 ;         
    end
    GR = result_KO.vBM ; Ex = sum(result_KO.Ex) ;
    result_array_Gr_KO_org_2(i) = GR(2)/GR_WT_org_2 ;
    result_array_Gr_WT_org_1(i) = GR(1)/GR_WT_org_1 ;
    result_array_Ex_org_2(i) = Ex ;
    result_array_abundance_org_2(i) = result_KO.BM(1) ;
    result_array_GRmax_KO_org_2(i) = result_KO.GRmax ; 
    
    catch
        fault_gene_ind_org_2 = [fault_gene_ind_org_2; targ_gene_ind] ; 
    end
end

gene_ind = single_lethal_genes_ind_org_2(find(result_array_org_2)) ; 
gr_KO_org_2 = result_array_Gr_KO_org_2(find(result_array_org_2)) ; 
gr_WT_org_1 = result_array_Gr_WT_org_1(find(result_array_org_2)) ; 
ex_flux = result_array_Ex_org_2(find(result_array_org_2)) ;
grmax = result_array_GRmax_KO_org_2(find(result_array_org_2)) ;
genes = model_2_WT.genes(gene_ind) ; 
abundance_KO = result_array_abundance_org_2(find(result_array_org_2)) ;
Table_org_2 = table(gene_ind,genes, (gr_KO_org_2'*GR_WT_org_2)', (gr_WT_org_1'*GR_WT_org_1)',gr_KO_org_2,grmax,abundance_KO,ex_flux, 'VariableNames',...
             {'Gene Index','Gene', 'KO Growth rate' , 'WT Growth rate', 'Ratio (gr_KO/gr_WT)','GRmax','KO Abundance','Sum of Ex Flux'}) ;
writetable(Table_org_2, excel_name,'Sheet',2) ; 
% Storing info into tables
result_struct_org_2 = table2struct(Table_org_2) ; 
final_results(1).Org_1_WT_Org2_KO_info_table = result_struct_org_2 ;
final_results(1).Org_1_WT_Org2_KO_info.Growth_rate_ratio_Org_1 = result_array_Gr_WT_org_1 ; 
final_results(1).Org_1_WT_Org2_KO_info.Growth_rate_ratio_Org_2 = result_array_Gr_KO_org_2 ; 
final_results(1).Org_1_WT_Org2_KO_info.Exchange_flux_sum = result_array_Ex_org_2 ; 
final_results(1).Org_1_WT_Org2_KO_info.Abundance_Org_1 = result_array_abundance_org_2 ; 


figure
bar(single_lethal_genes_ind_org_2, result_array_Gr_KO_org_2)
xlabel('Index of Gene in org_2 deleted') ; 
ylabel('Growth Rate ratio for Mutant (org_2) in Comm Model')
saveas(gcf, "glucose_org_2_KO_gr_ratio.png");

figure
bar(single_lethal_genes_ind_org_2, result_array_Gr_WT_org_1)
xlabel('Index of Gene in org_2 deleted') ; 
ylabel('Growth Rate ratio for WT (org_1) in Comm Model')
saveas(gcf, "glucose_model_1_WT_gr_ratio.png");

figure
bar(single_lethal_genes_ind_org_2, result_array_Ex_org_2)
xlabel('Index of Gene in org_2 deleted') ; 
ylabel('Exchange Fluxes (mean)')
saveas(gcf, "glucose_org_2_KO_Ex.png");

figure
bar(single_lethal_genes_ind_org_2, result_array_abundance_org_2)
xlabel('Index of Gene in org_2 deleted') ; 
ylabel('org_1 Abundance')
saveas(gcf, "glucose_org_2_KO_Abundance.png");


%% Checking performance of Communities with KO-KO combination, successful KOs from previous two simulations

target_genes_ind_org_1 = single_lethal_genes_ind_org_1(result_array_org_1 == 1) ;  % ind of genes to be KO in WT org_1
target_genes_ind_org_2 = single_lethal_genes_ind_org_2(result_array_org_2 == 1) ;  % ind of genes to be KO in WT org_2

% Creating arrays to store important results
num_succ_KO_org_1 = length(target_genes_ind_org_1) ; % Number of KO strains of org_1 that could grow
num_succ_KO_org_2 = length(target_genes_ind_org_2) ; % Number of KO strains of org_2 that could grow
total_combinations = num_succ_KO_org_1*num_succ_KO_org_2 ;      % total combinations possible

%auxotropic_pool_succ_KO = cell(total_combinations,3) ;          % array to store succesful auxotrophic strains
del_gene_array_KO_org_1_KO_org_2 = cell(num_succ_KO_org_1,num_succ_KO_org_2) ;  % array to store succesful gene deletions (gene id)
del_all_gene_KO_org_1 = ones(num_succ_KO_org_1,num_succ_KO_org_2)  ;                  % array to store all org_1 deleted gene info
del_all_gene_KO_org_2 = ones(num_succ_KO_org_1,num_succ_KO_org_2)  ;                  % array to store all org_2 deleted gene info

result_abundance_org_1 = zeros(total_combinations,1) ;    % Array to store org_1 abundances
result_abundance_org_2 = zeros(total_combinations,1) ;    % Array to store org_2 abundances
abundance_array = zeros(num_succ_KO_org_1,num_succ_KO_org_2) ; %  to store abundance in matrix

result_array_gr_KO_org_1_KO_org_2 = zeros(num_succ_KO_org_1,num_succ_KO_org_2) ; % ratio of growth rate of KO org_1 with KO org_2 to WT org_1 
result_array_gr_KO_org_2_KO_org_1 = zeros(num_succ_KO_org_1,num_succ_KO_org_2) ; % ratio of growth rate of KO org_2 with KO org_1 to WT org_2
result_array_gr_KO_comm = zeros(num_succ_KO_org_1,num_succ_KO_org_2) ; % GRmax of community 

result_array_Ex_KO_org_1_KO_org_2 = zeros(num_succ_KO_org_1,num_succ_KO_org_2) ;  % Export fluxes
   
error_comb_array = cell(num_succ_KO_org_1,num_succ_KO_org_2)   ;

%count = 1 ; 

parfor m = 1:num_succ_KO_org_1
    
    for j = 1:num_succ_KO_org_2
        restoreEnvironment(environment);
        disp([m,j])
%         if ismember(m,[10,11,12,24,25,26,27,28,30,34,35,36,41,47, ...
%                 49,53,54,59,60]) % &&  ismember(j,[1,11,28])%
%                 ismember(m,[16,17,23]) &&  ismember(j,[71,81,108]) % %
%                 skipping any indices means that I encountered some issue with the solver and it used to run infinitely so skipped them for now.  
%               %ind = ind + 1 ;
%               continue
%          end
%          if ismember(m,[10,11,12,24])  &&  ismember(j,[1,11,16,17,28,44,46,52,53])% ismember(m,[16,17,23]) &&  ismember(j,[71,81,108]) %
%               %ind = ind + 1 ;
%               continue
%          end
         
         
        tar_gene_ind_org_1 = target_genes_ind_org_1(m) ; 
        tar_gene_ind_org_2 = target_genes_ind_org_2(j) ; 
        
        ko_gene_org_1 = strcat(model_1_WT.genes(tar_gene_ind_org_1) , {'_model_1'}) ; 
        ko_gene_org_2 = strcat(model_2_WT.genes(tar_gene_ind_org_2) , {'_model_2'}) ; 
        ko_comm_model = deleteModelGenes(comm_model, [ko_gene_org_1,ko_gene_org_2]) ; 
        
        try 
        [sol_KO, result_KO] = SteadyCom(ko_comm_model) ; % Applying Steady Come to KO community model
        
         if result_KO.stat == "optimal"
            if result_KO.vBM(1) > 1e-5 && result_KO.vBM(2) > 1e-5
                %auxotropic_pool_succ_KO{ind} = {comm_models{1}, org_1_KO, org_2_KO} ;% causes memory issue
                del_gene_array_KO_org_1_KO_org_2{m,j} = {tar_gene_ind_org_1,tar_gene_ind_org_2} ;                
            end  
        end
        
        result_array_gr_KO_org_1_KO_org_2(m,j) = result_KO.vBM(1)/GR_WT_org_1 ; 
        result_array_gr_KO_org_2_KO_org_1(m,j) = result_KO.vBM(2)/GR_WT_org_2 ; 
        result_array_gr_KO_comm(m,j) = result_KO.GRmax ; 
        result_array_Ex_KO_org_1_KO_org_2(m,j) = sum(result_KO.Ex) ; 
        del_all_gene_KO_org_1(m,j) = tar_gene_ind_org_1 ;  
        del_all_gene_KO_org_2(m,j) = tar_gene_ind_org_2  ;       
        abundance_array(m,j) = result_KO.BM(1)  ;        
        
        catch 
            error_comb_array{m,j} = {tar_gene_ind_org_1,tar_gene_ind_org_2} ;  
            del_all_gene_KO_org_1(m,j) = tar_gene_ind_org_1 ;  
            del_all_gene_KO_org_2(m,j) = tar_gene_ind_org_2  ;
             
        end
        
    end
end

figure
imagesc(result_array_gr_KO_org_1_KO_org_2)
colorbar
xlabel('Gene Ind of deleted org_2 Gene in target gene array') ; 
ylabel('Gene Ind of deleted org_1 Gene in target gene array') ; 
title([['Growth rate ratio for KO org_1 w.r.t WT org_1'],['In community of KO strains']]) ;
saveas(gcf, "glucose_comm_gr_ratio_org_1.png");

figure
imagesc(result_array_gr_KO_org_2_KO_org_1)
colorbar
xlabel('Gene Ind of deleted org_2 Gene in target gene array') ; 
ylabel('Gene Ind of deleted org_1 Gene in target gene array') ; 
title([['Growth rate ratio for KO org_2 w.r.t WT org_2'],['In community of KO strains']]) ;
saveas(gcf, "glucose_comm_gr_ratio_org_2.png");

figure
imagesc(result_array_gr_KO_comm)
colorbar
xlabel('Gene Ind of deleted org_2 Gene in target gene array') ; 
ylabel('Gene Ind of deleted org_1 Gene in target gene array') ; 
title([['Max Growth rate (GRmax)'],['In community of KO strains']]) ;
saveas(gcf, "glucose_comm_GRmax.png");

        
figure
imagesc(result_array_Ex_KO_org_1_KO_org_2)
colorbar
xlabel('Gene Ind of deleted org_2 Gene in target gene array') ; 
ylabel('Gene Ind of deleted org_1 Gene in target gene array') ; 
title([['Total Ex Flux'],['In community of KO strains']]) ;
saveas(gcf, "glucose_comm_Ex.png");

figure
imagesc(abundance_array)
colorbar
xlabel('Gene Ind of deleted org_2 Gene in target gene array') ; 
ylabel('Gene Ind of deleted org_1 Gene in target gene array') ; 
title([['Abundance (org_1)'],['In community of KO strains']]) ;
saveas(gcf, "glucose_comm_abundance.png");

result_array_gr_KO_org_1_KO_org_2_trans = result_array_gr_KO_org_1_KO_org_2' ; 
result_ratio_gr_org_1 = result_array_gr_KO_org_1_KO_org_2_trans(:) ; 

result_array_gr_KO_org_2_KO_org_1_trans = result_array_gr_KO_org_2_KO_org_1' ; 
result_ratio_gr_org_2 = result_array_gr_KO_org_2_KO_org_1_trans(:) ; 

result_array_gr_KO_comm_trans = result_array_gr_KO_comm' ; 
result_grmax = result_array_gr_KO_comm_trans(:) ; 

result_array_Ex_KO_org_1_KO_org_2_trans = result_array_Ex_KO_org_1_KO_org_2' ; 
result_ex_flux = result_array_Ex_KO_org_1_KO_org_2_trans(:) ; 

 

del_all_gene_KO_org_1_trans = del_all_gene_KO_org_1' ; 
del_all_gene_KO_org_1 = del_all_gene_KO_org_1_trans(:) ; 

del_all_gene_KO_org_2_trans = del_all_gene_KO_org_2' ; 
del_all_gene_KO_org_2 = del_all_gene_KO_org_2_trans(:) ; 


del_gene_org_1_ind = model_1_WT.genes(del_all_gene_KO_org_1) ; 
del_gene_org_2_ind = model_2_WT.genes(del_all_gene_KO_org_2) ; 


ind = 1 : total_combinations ; 

result_abundance_org_1_trans = abundance_array' ; 
result_abundance = result_abundance_org_1_trans(:) ; 

Table_succ = table(ind', del_gene_org_1_ind, del_gene_org_2_ind, result_ratio_gr_org_1, result_ratio_gr_org_2, ...
    result_grmax, result_ex_flux, result_abundance,'VariableNames',...
    {'Sr.no', 'org_1 del gene','org_2 del gene','org_1 Ratio (gr_KO/gr_WT)', 'org_2 Ratio (gr_KO/gr_WT)', 'GRmax', 'Ex Flux', 'org_1 abundace'}); 
writetable(Table_succ, excel_name,'Sheet',3) ;

% Storing results to final structure 
result_struct_succ_KO = table2struct(Table_succ) ; 
final_results(1).Successful_KOs_info_table = result_struct_succ_KO ;
final_results(1).Successful_KOs_info.Growth_rate_ratio_Org_1 = result_array_gr_KO_org_1_KO_org_2 ; 
final_results(1).Successful_KOs_info.Growth_rate_ratio_Org_2 = result_array_gr_KO_org_2_KO_org_1 ; 
final_results(1).Successful_KOs_info.Max_Growth_rate = result_array_gr_KO_comm ; 
final_results(1).Successful_KOs_info.Exchange_Fluxes = result_array_Ex_KO_org_1_KO_org_2 ; 
final_results(1).Successful_KOs_info.Abundances_Org_1 = abundance_array ;


table_ab = table(['org_1' ; 'org_2'],wild_type_abundances,'VariableNames',{'Org', 'WT-WT Abundances'}) ;
writetable(table_ab, excel_name,'Sheet',4) ;
toc
end