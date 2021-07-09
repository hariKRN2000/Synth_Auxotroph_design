close all ;clear; clc; 
tic
% Setting up solver
%changeCobraSolver('gurobi')


% Loading Models
ecoli_WT = readCbModel('iML1515.mat');
yeast_WT = readCbModel('iMM904.mat') ;

% Setting Anaerobic Condition : 


all_models = {ecoli_WT,yeast_WT} ; % Creating model list for WT

[comm_models,pairedModelInfo] = create_community(all_models) ; % Creating community model for WT 

%comm_models{1} = changeRxnBounds(comm_models{1},'EX_o2_e(u)',-0,'l'); % Set maximum oxygen uptake
%comm_models{1} = changeRxnBounds(comm_models{1},'EX_glc__D_e(u)',0,'l'); %To remove glucose from medium
%comm_models{1} = changeRxnBounds(comm_models{1},'EX_malt_e(u)',-10,'l'); %To change carbon source

[sol_WT, result_WT] = SteadyCom(comm_models{1}) ; % Applying SteadyCom to WT community model 
GR_WT = result_WT.GRmax ; % Growth rate for WT 
GR_WT_ecoli = result_WT.vBM(1)  % Growth Rate for WT Ecoli
GR_WT_yeast = result_WT.vBM(2)  % Growth Rate for WT Yeast
wild_type_abundances = result_WT.BM % Abundances

% Finding Lethal genes for WT Yeast
[grRatio_single_yeast, grRateKO_single_yeast] = singleGeneDeletion(yeast_WT);                               
single_lethal_genes_ind_yeast = find(grRatio_single_yeast < 0.05 ) ;
% Finding Lethal genes for WT Ecoli
[grRatio_single_ecoli, grRateKO_single_ecoli] = singleGeneDeletion(ecoli_WT);                               
single_lethal_genes_ind_ecoli = find(grRatio_single_ecoli < 0.05 ) ;

% Creating arrays to store important results
result_array_ecoli = zeros(length(single_lethal_genes_ind_ecoli),1) ; % array noting which KO grows with WT 
result_array_yeast = zeros(length(single_lethal_genes_ind_yeast),1) ; % array noting which KO grows with WT 

%auxotropic_pool_ecoli = cell(length(single_lethal_genes_ind_ecoli),3) ; % array to store succesful auxotrophic strains
%auxotropic_pool_yeast = cell(length(single_lethal_genes_ind_yeast),3) ; % array to store succesful auxotrophic strains

result_array_Gr_KO_ecoli = zeros(length(single_lethal_genes_ind_ecoli),1) ; % ratio of growth rate of KO to WT 
result_array_Gr_KO_yeast = zeros(length(single_lethal_genes_ind_yeast),1) ; % ratio of growth rate of KO to WT 

result_array_Gr_WT_ecoli = zeros(length(single_lethal_genes_ind_yeast),1) ; % ratio of growth rate of WT with KO to WT 
result_array_Gr_WT_yeast = zeros(length(single_lethal_genes_ind_ecoli),1) ; % ratio of growth rate of WT with KO to WT 

result_array_Ex_ecoli = zeros(length(single_lethal_genes_ind_ecoli),1) ; % Export fluxes where KO is ecoli
result_array_Ex_yeast = zeros(length(single_lethal_genes_ind_yeast),1) ; % Export fluxes where KO is yeast

result_array_abundance_ecoli = zeros(length(single_lethal_genes_ind_ecoli),1) ; % Abundance of ecoli where KO is ecoli
result_array_abundance_yeast = zeros(length(single_lethal_genes_ind_yeast),1) ; % Abundance of ecoli where KO is yeast

result_array_GRmax_KO_ecoli = zeros(length(single_lethal_genes_ind_ecoli),1) ;
result_array_GRmax_KO_yeast = zeros(length(single_lethal_genes_ind_yeast),1) ;
fault_gene_ind_ecoli = [] ;
fault_gene_ind_yeast = [] ;
 

% Simulating models for KO Ecoli : 
parfor i = 1:length(single_lethal_genes_ind_ecoli) 
    disp(single_lethal_genes_ind_ecoli(i))
    disp(i)
    
%     if ismember(i,[41]) % &&  ismember(j,error_ind_gene_yeast)
%              %ind = ind + 1 ;
%              continue
%     end
    
%     model_1 = deleteModelGenes(ecoli_WT, ecoli_WT.genes(single_lethal_genes_ind_ecoli(i)));
%     all_models = {model_1 , yeast_WT} ;
%     [comm_models,pairedModelInfo] = create_community(all_models) ; % Creating community model for WT
    targ_gene_ind = single_lethal_genes_ind_ecoli(i) ; 

    ko_gene_ecoli = strcat(ecoli_WT.genes(targ_gene_ind) , {'_model_1'}) ; 
  
    ko_comm_model = deleteModelGenes(comm_models{1}, ko_gene_ecoli) ;
    %ko_comm_model = changeRxnBounds(ko_comm_model,'EX_o2_e(u)',-0,'l'); % Set maximum oxygen uptake
    %comm_models{1} = changeRxnBounds(comm_models{1},'EX_o2_e(u)',-0,'l'); % Set maximum oxygen uptake
    
    try
    [sol_KO, result_KO] = SteadyCom(ko_comm_model) ; % Applying Steady Come to WT community model
    
    if result_KO.vBM(1) > 0 
        %auxotropic_pool_ecoli{i} = {comm_models{1}, model_1, result_KO} ;
        result_array_ecoli(i) = 1 ;         
    end
    GR = result_KO.vBM ; Ex = sum(result_KO.Ex) ;
    result_array_Gr_KO_ecoli(i) = GR(1)/GR_WT_ecoli ;
    result_array_Gr_WT_yeast(i) = GR(2)/GR_WT_yeast ;
    result_array_Ex_ecoli(i) = Ex ;
    result_array_abundance_ecoli(i) = result_KO.BM(1) ; 
    result_array_GRmax_KO_ecoli(i) = result_KO.GRmax ;
    
    catch
        fault_gene_ind_ecoli = [fault_gene_ind_ecoli; targ_gene_ind] ;
    end
    
end

gene_ind = single_lethal_genes_ind_ecoli(find(result_array_ecoli)) ; 
gr_KO_ecoli = result_array_Gr_KO_ecoli(find(result_array_ecoli)) ; 
gr_WT_yeast = result_array_Gr_WT_yeast(find(result_array_ecoli)) ; 
ex_flux = result_array_Ex_ecoli(find(result_array_ecoli)) ;
grmax = result_array_GRmax_KO_ecoli(find(result_array_ecoli)) ;
genes = ecoli_WT.genes(gene_ind) ; 
abundance_KO = result_array_abundance_ecoli(find(result_array_ecoli)) ;
Table_ecoli = table(gene_ind,genes, gr_KO_ecoli'*GR_WT_ecoli, gr_WT_yeast'*GR_WT_yeast,gr_KO_ecoli',grmax,abundance_KO,ex_flux, 'VariableNames',...
             {'Gene Index','Gene', 'KO Growth rate' , 'WT Growth rate', 'Ratio (gr_KO/gr_WT)','GRmax','KO Abundance','Sum of Ex Flux'}) ;
writetable(Table_ecoli, 'test_KO_glucose.xls','Sheet',1) ; 

figure
bar(single_lethal_genes_ind_ecoli, result_array_Gr_KO_ecoli)
xlabel('Index of Gene in Ecoli deleted') ; 
ylabel('Growth Rate ratio for Mutant (Ecoli) in Comm Model')
saveas(gcf, "glucose_ecoli_KO_gr_ratio.png");

figure
bar(single_lethal_genes_ind_ecoli, result_array_Gr_WT_yeast)
xlabel('Index of Gene in Ecoli deleted') ; 
ylabel('Growth Rate ratio for WT (Yeast) in Comm Model')
saveas(gcf, "glucose_yeast_WT_gr_ratio.png");


figure
bar(single_lethal_genes_ind_ecoli, result_array_Ex_ecoli)
xlabel('Index of Gene in Ecoli deleted') ; 
ylabel('Exchange Fluxes (mean)')
saveas(gcf, "glucose_ecoli_KO_Ex.png");

figure
bar(single_lethal_genes_ind_ecoli, result_array_abundance_ecoli)
xlabel('Index of Gene in Ecoli deleted') ; 
ylabel('Ecoli Abundance')
saveas(gcf, "glucose_ecoli_KO_Abundance.png");


% Simulating models for KO Yeast : 
parfor i = 1:length(single_lethal_genes_ind_yeast) 
    
%     model_1 = deleteModelGenes(yeast_WT, yeast_WT.genes(single_lethal_genes_ind_yeast(i)));
%     all_models = {model_1 , ecoli_WT} ;
%     [comm_models,pairedModelInfo] = create_community(all_models) ; % Creating community model for WT
    
    targ_gene_ind = single_lethal_genes_ind_yeast(i) ;

    ko_gene_yeast = strcat(yeast_WT.genes(targ_gene_ind) , {'_model_2'}) ; 
  
    ko_comm_model = deleteModelGenes(comm_models{1}, ko_gene_yeast) ;
    
    %ko_comm_model = changeRxnBounds(ko_comm_model,'EX_o2_e(u)',-0,'l'); % Set maximum oxygen uptake
    %comm_models{1} = changeRxnBounds(comm_models{1},'EX_o2_e(u)',-0,'l'); % Set maximum oxygen uptake
   
    try
    [sol_KO, result_KO] = SteadyCom(ko_comm_model) ; % Applying Steady Come to WT community model

    
    if result_KO.vBM(1) > 0 
        %auxotropic_pool_yeast{i} = {comm_models{1}, model_1, result_KO} ;
        result_array_yeast(i) = 1 ;         
    end
    GR = result_KO.vBM ; Ex = sum(result_KO.Ex) ;
    result_array_Gr_KO_yeast(i) = GR(2)/GR_WT_yeast ;
    result_array_Gr_WT_ecoli(i) = GR(1)/GR_WT_ecoli ;
    result_array_Ex_yeast(i) = Ex ;
    result_array_abundance_yeast(i) = result_KO.BM(1) ;
    result_array_GRmax_KO_yeast(i) = result_KO.GRmax ; 
    
    catch
        fault_gene_ind_yeast = [fault_gene_ind_yeast; targ_gene_ind] ; 
    end
end

gene_ind = single_lethal_genes_ind_yeast(find(result_array_yeast)) ; 
gr_KO_yeast = result_array_Gr_KO_yeast(find(result_array_yeast)) ; 
gr_WT_ecoli = result_array_Gr_WT_ecoli(find(result_array_yeast)) ; 
ex_flux = result_array_Ex_yeast(find(result_array_yeast)) ;
grmax = result_array_GRmax_KO_yeast(find(result_array_yeast)) ;
genes = yeast_WT.genes(gene_ind) ; 
abundance_KO = result_array_abundance_yeast(find(result_array_yeast)) ;
Table_yeast = table(gene_ind,genes, gr_KO_yeast'*GR_WT_yeast, gr_WT_ecoli'*GR_WT_ecoli,gr_KO_yeast',grmax,abundance_KO,ex_flux, 'VariableNames',...
             {'Gene Index','Gene', 'KO Growth rate' , 'WT Growth rate', 'Ratio (gr_KO/gr_WT)','GRmax','KO Abundance','Sum of Ex Flux'}) ;
writetable(Table_yeast, 'test_KO_glucose.xls','Sheet',2) ; 

figure
bar(single_lethal_genes_ind_yeast, result_array_Gr_KO_yeast)
xlabel('Index of Gene in Yeast deleted') ; 
ylabel('Growth Rate ratio for Mutant (Yeast) in Comm Model')
saveas(gcf, "glucose_yeast_KO_gr_ratio.png");

figure
bar(single_lethal_genes_ind_yeast, result_array_Gr_WT_ecoli)
xlabel('Index of Gene in yeast deleted') ; 
ylabel('Growth Rate ratio for WT (Ecoli) in Comm Model')
saveas(gcf, "glucose_ecoli_WT_gr_ratio.png");

figure
bar(single_lethal_genes_ind_yeast, result_array_Ex_yeast)
xlabel('Index of Gene in Yeast deleted') ; 
ylabel('Exchange Fluxes (mean)')
saveas(gcf, "glucose_yeast_KO_Ex.png");

figure
bar(single_lethal_genes_ind_yeast, result_array_abundance_yeast)
xlabel('Index of Gene in Yeast deleted') ; 
ylabel('Ecoli Abundance')
saveas(gcf, "glucose_yeast_KO_Abundance.png");








% Testing succesfull KOs



target_genes_ind_ecoli = single_lethal_genes_ind_ecoli(result_array_ecoli == 1) ;  % ind of genes to be KO in WT Ecoli
target_genes_ind_yeast = single_lethal_genes_ind_yeast(result_array_yeast == 1) ;  % ind of genes to be KO in WT Yeast

% Creating arrays to store important results
num_succ_KO_ecoli = length(target_genes_ind_ecoli) ; % Number of KO strains of ecoli that could grow
num_succ_KO_yeast = length(target_genes_ind_yeast) ; % Number of KO strains of yeast that could grow
total_combinations = num_succ_KO_ecoli*num_succ_KO_yeast ;      % total combinations possible

%auxotropic_pool_succ_KO = cell(total_combinations,3) ;          % array to store succesful auxotrophic strains
del_gene_array_KO_ecoli_KO_yeast = zeros(total_combinations,2) ;  % array to store succesful gene deletions (gene id)
del_all_gene_KO = ones(total_combinations,2) ;                  % array to store all deleted gene info

result_abundance_ecoli = zeros(total_combinations,1) ;    % Array to store ecoli abundances
result_abundance_yeast = zeros(total_combinations,1) ;    % Array to store yeast abundances
abundance_array = zeros(num_succ_KO_ecoli,num_succ_KO_yeast) ; %  to store abundance in matrix

result_array_gr_KO_ecoli_KO_yeast = zeros(num_succ_KO_ecoli,num_succ_KO_yeast) ; % ratio of growth rate of KO ecoli with KO yeast to WT ecoli 
result_array_gr_KO_yeast_KO_ecoli = zeros(num_succ_KO_ecoli,num_succ_KO_yeast) ; % ratio of growth rate of KO yeast with KO ecoli to WT yeast
result_array_gr_KO_comm = zeros(num_succ_KO_ecoli,num_succ_KO_yeast) ; % GRmax of community 

result_array_Ex_KO_ecoli_KO_yeast = zeros(num_succ_KO_ecoli,num_succ_KO_yeast) ;  % Export fluxes

%  error_ind_gene_ecoli = [2,8,19,28,44,57,58] ; 
%  error_ind_gene_yeast = [102,104,103,105,109,110] ; 
%  
%  error_ind_gene_ecoli = [86] ; 
%  error_ind_gene_yeast = [102,105,109,110] ; 
%  
error_comb_array = zeros(total_combinations,2) ;

count = 1 ; 
parfor i = 1:num_succ_KO_ecoli
    parfor j = 1:num_succ_KO_yeast
        disp([i,j])
%         if ismember(i,error_ind_gene_ecoli) % &&  ismember(j,error_ind_gene_yeast)
%              ind = ind + 1 ;
%              continue
%         end
        tar_gene_ind_ecoli = target_genes_ind_ecoli(i) ; 
        tar_gene_ind_yeast = target_genes_ind_yeast(j) ; 
        %ecoli_KO = deleteModelGenes(ecoli_WT, ecoli_WT.genes(tar_gene_ind_ecoli));
        %yeast_KO = deleteModelGenes(yeast_WT, yeast_WT.genes(tar_gene_ind_yeast));
        %all_models = {ecoli_KO, yeast_KO} ;
        %[comm_models,pairedModelInfo] = create_community(all_models) ; % Creating community model
        ko_gene_ecoli = strcat(ecoli_WT.genes(tar_gene_ind_ecoli) , {'_model_1'}) ; 
        ko_gene_yeast = strcat(yeast_WT.genes(tar_gene_ind_yeast) , {'_model_2'}) ; 
        ko_comm_model = deleteModelGenes(comm_models{1}, [ko_gene_ecoli,ko_gene_yeast]) ; 
        %ko_comm_model = changeRxnBounds(ko_comm_model,'EX_o2_e(u)',-0,'l'); % Set maximum oxygen uptake
        
        %fbasol = optimizeCbModel(ko_comm_model) ;% First Testing for FBA
       
        try 
        [sol_KO, result_KO] = SteadyCom(ko_comm_model) ; % Applying Steady Come to KO community model
        
         if result_KO.stat == "optimal"
            if result_KO.vBM(1) > 1e-5 && result_KO.vBM(2) > 1e-5
                %auxotropic_pool_succ_KO{ind} = {comm_models{1}, ecoli_KO, yeast_KO} ;% causes memory issue
                del_gene_array_KO_ecoli_KO_yeast(count,1) = tar_gene_ind_ecoli ;
                del_gene_array_KO_ecoli_KO_yeast(count,2) = tar_gene_ind_yeast ;
            end  
        end
        
        result_array_gr_KO_ecoli_KO_yeast(i,j) = result_KO.vBM(1)/GR_WT_ecoli ; 
        result_array_gr_KO_yeast_KO_ecoli(i,j) = result_KO.vBM(2)/GR_WT_yeast ; 
        result_array_gr_KO_comm(i,j) = result_KO.GRmax ; 
        result_array_Ex_KO_ecoli_KO_yeast(i,j) = sum(result_KO.Ex) ; 
        del_all_gene_KO(count,1) = tar_gene_ind_ecoli ; 
        del_all_gene_KO(count,2) = tar_gene_ind_yeast ; 
         
        result_abundance_ecoli(count) = result_KO.BM(1)  ; 
        result_abundance_yeast(count) = result_KO.BM(2)  ; 
        abundance_array(i,j) = result_KO.BM(1)  ;
        %ind = ind + 1 ;
        
        catch 
            error_comb_array(count,1) = tar_gene_ind_ecoli ; 
            error_comb_array(count,2) = tar_gene_ind_yeast ;
            del_all_gene_KO(count,1) = tar_gene_ind_ecoli ; 
            del_all_gene_KO(count,2) = tar_gene_ind_yeast ;            
            %ind = ind + 1 ; 
        end
        count = count + 1 ;
    end
end

figure
imagesc(result_array_gr_KO_ecoli_KO_yeast)
colorbar
xlabel('Gene Ind of deleted Yeast Gene in target gene array') ; 
ylabel('Gene Ind of deleted Ecoli Gene in target gene array') ; 
title([['Growth rate ratio for KO Ecoli w.r.t WT Ecoli'],['In community of KO strains']]) ;
saveas(gcf, "glucose_comm_gr_ratio_ecoli.png");

figure
imagesc(result_array_gr_KO_yeast_KO_ecoli)
colorbar
xlabel('Gene Ind of deleted Yeast Gene in target gene array') ; 
ylabel('Gene Ind of deleted Ecoli Gene in target gene array') ; 
title([['Growth rate ratio for KO Yeast w.r.t WT Yeast'],['In community of KO strains']]) ;
saveas(gcf, "glucose_comm_gr_ratio_yeast.png");

figure
imagesc(result_array_gr_KO_comm)
colorbar
xlabel('Gene Ind of deleted Yeast Gene in target gene array') ; 
ylabel('Gene Ind of deleted Ecoli Gene in target gene array') ; 
title([['Max Growth rate (GRmax)'],['In community of KO strains']]) ;
saveas(gcf, "glucose_comm_GRmax.png");

        
figure
imagesc(result_array_Ex_KO_ecoli_KO_yeast)
colorbar
xlabel('Gene Ind of deleted Yeast Gene in target gene array') ; 
ylabel('Gene Ind of deleted Ecoli Gene in target gene array') ; 
title([['Total Ex Flux'],['In community of KO strains']]) ;
saveas(gcf, "glucose_comm_Ex.png");

figure
imagesc(abundance_array)
colorbar
xlabel('Gene Ind of deleted Yeast Gene in target gene array') ; 
ylabel('Gene Ind of deleted Ecoli Gene in target gene array') ; 
title([['Abundance (Ecoli)'],['In community of KO strains']]) ;
saveas(gcf, "glucose_comm_abundance.png");

result_array_gr_KO_ecoli_KO_yeast_trans = result_array_gr_KO_ecoli_KO_yeast' ; 
result_ratio_gr_ecoli = result_array_gr_KO_ecoli_KO_yeast_trans(:) ; 

result_array_gr_KO_yeast_KO_ecoli_trans = result_array_gr_KO_yeast_KO_ecoli' ; 
result_ratio_gr_yeast = result_array_gr_KO_yeast_KO_ecoli_trans(:) ; 

result_array_gr_KO_comm_trans = result_array_gr_KO_comm' ; 
result_grmax = result_array_gr_KO_comm_trans(:) ; 

result_array_Ex_KO_ecoli_KO_yeast_trans = result_array_Ex_KO_ecoli_KO_yeast' ; 
result_ex_flux = result_array_Ex_KO_ecoli_KO_yeast_trans(:) ; 

del_gene_ecoli_ind = ecoli_WT.genes(del_all_gene_KO(:,1)) ; 
del_gene_yeast_ind = yeast_WT.genes(del_all_gene_KO(:,2)) ; 

ind = 1 : total_combinations ; 

% result_abundance_ecoli_trans = result_abundance_ecoli' ; 
% result_abundance = result_abundance_ecoli_trans(:) ; 

Table_succ = table(ind', del_gene_ecoli_ind, del_gene_yeast_ind, result_ratio_gr_ecoli, result_ratio_gr_yeast, ...
    result_grmax, result_ex_flux, result_abundance_ecoli,'VariableNames',...
    {'Sr.no', 'Ecoli del gene','Yeast del gene','Ecoli Ratio (gr_KO/gr_WT)', 'Yeast Ratio (gr_KO/gr_WT)', 'GRmax', 'Ex Flux', 'Ecoli abundace'}); 
writetable(Table_succ, 'test_KO_glucose.xls','Sheet',3) ;

table_ab = table(['Ecoli' ; 'Yeast'],wild_type_abundances,'VariableNames',{'Org', 'WT-WT Abundances'}) ;
writetable(table_ab, 'test_KO_glucose.xls','Sheet',4) ;
toc