clear; clc; 
% Loading models 
ecoli_WT = readCbModel('iML1515.mat');
yeast_WT = readCbModel('iMM904.mat') ;

% Storing WT community model
all_models = {ecoli_WT,yeast_WT} ; % Creating model list for WT

[comm_models,pairedModelInfo] = create_community(all_models) ; % Creating community model for WT 

[sol_WT, result_WT] = SteadyCom(comm_models{1}) ; % Applying SteadyCom to WT community model 
GR_WT = result_WT.GRmax ; % Growth rate for WT 
GR_WT_ecoli = result_WT.vBM(1) ; % Growth Rate for WT Ecoli
GR_WT_yeast = result_WT.vBM(2) ; % Growth Rate for WT Yeast

% Importing necessary data from saved tables
Table_yeast = readtable('iMM904_iML1515_result_table.csv') ; 
Table_ecoli = readtable('iML1515_iMM904_result_table.csv') ;

% Creating arrays


target_genes_ind_ecoli = Table_ecoli{:,1} ;  % ind of genes to be KO in WT Ecoli
target_genes_ind_yeast = Table_yeast{:,1} ;  % ind of genes to be KO in WT Yeast

% Creating arrays to store important results
num_succ_KO_ecoli = length(target_genes_ind_ecoli) ; % Number of KO strains of ecoli that could grow
num_succ_KO_yeast = length(target_genes_ind_yeast) ; % Number of KO strains of yeast that could grow
total_combinations = num_succ_KO_ecoli*num_succ_KO_yeast ;      % total combinations possible

%auxotropic_pool_succ_KO = cell(total_combinations,3) ;          % array to store succesful auxotrophic strains
del_gene_array_KO_ecoli_KO_yeast = zeros(total_combinations,2) ;  % array to store succesful gene deletions (gene id)
del_all_gene_KO = zeros(total_combinations,2) ;                  % array to store all deleted gene info

result_abundance_ecoli = zeros(total_combinations,1) ;


result_array_gr_KO_ecoli_KO_yeast = zeros(num_succ_KO_ecoli,num_succ_KO_yeast) ; % ratio of growth rate of KO ecoli with KO yeast to WT ecoli 
result_array_gr_KO_yeast_KO_ecoli = zeros(num_succ_KO_ecoli,num_succ_KO_yeast) ; % ratio of growth rate of KO yeast with KO ecoli to WT yeast
result_array_gr_KO_comm = zeros(num_succ_KO_ecoli,num_succ_KO_yeast) ; % GRmax of community 

result_array_Ex_KO_ecoli_KO_yeast = zeros(num_succ_KO_ecoli,num_succ_KO_yeast) ;  % Export fluxes

ind = 1 ; 
for i = 1:num_succ_KO_ecoli
    for j = 1:num_succ_KO_yeast
        tar_gene_ind_ecoli = target_genes_ind_ecoli(i) ; 
        tar_gene_ind_yeast = target_genes_ind_yeast(j) ; 
        %ecoli_KO = deleteModelGenes(ecoli_WT, ecoli_WT.genes(tar_gene_ind_ecoli));
        %yeast_KO = deleteModelGenes(yeast_WT, yeast_WT.genes(tar_gene_ind_yeast));
        %all_models = {ecoli_KO, yeast_KO} ;
        %[comm_models,pairedModelInfo] = create_community(all_models) ; % Creating community model
        ko_gene_ecoli = strcat(ecoli_WT.genes(tar_gene_ind_ecoli) , {'_model_1'}) ; 
        ko_gene_yeast = strcat(yeast_WT.genes(tar_gene_ind_yeast) , {'_model_2'}) ; 
        ko_comm_model = deleteModelGenes(comm_models{1}, [ko_gene_ecoli,ko_gene_yeast]) ; 
        
        [sol_KO, result_KO] = SteadyCom(ko_comm_model) ; % Applying Steady Come to KO community model
        
         if result_KO.stat == "optimal"
            if result_KO.vBM(1) > 1e-5 && result_KO.vBM(2) > 1e-5
                %auxotropic_pool_succ_KO{ind} = {comm_models{1}, ecoli_KO, yeast_KO} ;% causes memory issue
                del_gene_array_KO_ecoli_KO_yeast(ind,1) = tar_gene_ind_ecoli ;
                del_gene_array_KO_ecoli_KO_yeast(ind,2) = tar_gene_ind_yeast ;
            end  
        end
        
        result_array_gr_KO_ecoli_KO_yeast(i,j) = result_KO.vBM(1)/GR_WT_ecoli ; 
        result_array_gr_KO_yeast_KO_ecoli(i,j) = result_KO.vBM(2)/GR_WT_yeast ; 
        result_array_gr_KO_comm(i,j) = result_KO.GRmax ; 
        result_array_Ex_KO_ecoli_KO_yeast(i,j) = sum(result_KO.Ex) ; 
        del_all_gene_KO(ind,1) = tar_gene_ind_ecoli ; 
        del_all_gene_KO(ind,2) = tar_gene_ind_yeast ; 
        result_abundance_ecoli(ind) = result_KO.BM(1) ; 
        ind = ind + 1 ;
    end
end

figure
imagesc(result_array_gr_KO_ecoli_KO_yeast)
colorbar
xlabel('Gene Ind of deleted Yeast Gene in target gene array') ; 
ylabel('Gene Ind of deleted Ecoli Gene in target gene array') ; 
title([['Growth rate ratio for KO Ecoli w.r.t WT Ecoli'],['In community of KO strains']]) ;

figure
imagesc(result_array_gr_KO_yeast_KO_ecoli)
colorbar
xlabel('Gene Ind of deleted Yeast Gene in target gene array') ; 
ylabel('Gene Ind of deleted Ecoli Gene in target gene array') ; 
title([['Growth rate ratio for KO Yeast w.r.t WT Yeast'],['In community of KO strains']]) ;

figure
imagesc(result_array_gr_KO_comm)
colorbar
xlabel('Gene Ind of deleted Yeast Gene in target gene array') ; 
ylabel('Gene Ind of deleted Ecoli Gene in target gene array') ; 
title([['Max Growth rate (GRmax)'],['In community of KO strains']]) ;

        
figure
imagesc(result_array_Ex_KO_ecoli_KO_yeast)
colorbar
xlabel('Gene Ind of deleted Yeast Gene in target gene array') ; 
ylabel('Gene Ind of deleted Ecoli Gene in target gene array') ; 
title([['Total Ex Flux'],['In community of KO strains']]) ;

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

result_abundance_ecoli_trans = result_abundance_ecoli' ; 
result_abundance = result_abundance_ecoli_trans(:) ; 

Table = table(ind', del_gene_ecoli_ind, del_gene_yeast_ind, result_ratio_gr_ecoli, result_ratio_gr_yeast, ...
    result_grmax, result_ex_flux, result_abundance,'VariableNames',...
    {'Sr.no', 'Ecoli del gene','Yeast del gene','Ecoli Ratio (gr_KO/gr_WT)', 'Yeast Ratio (gr_KO/gr_WT)', 'GRmax', 'Ex Flux', 'Ecoli abundace'}) 

