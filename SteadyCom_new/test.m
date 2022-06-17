clear; clc; 
ecoli_WT = readCbModel('iML1515.mat');
yeast_WT = readCbModel('iMM904.mat') ;

all_models = {ecoli_WT,yeast_WT} ; % Creating model list for WT

[comm_models,pairedModelInfo] = create_community(all_models) ; % Creating community model for WT 
comm_models{1} = changeRxnBounds(comm_models{1},'EX_o2_e(u)',-0,'l'); % Set maximum oxygen uptake

[sol_WT, result_WT] = SteadyCom(comm_models{1}) ; % Applying Steady Come to WT community model 
GR_WT = result_WT.GRmax ; % Growth rate for WT 
GR_WT_ecoli = result_WT.vBM(1) ; % Growth Rate for WT Ecoli
GR_WT_yeast = result_WT.vBM(2) ; % Growth Rate for WT Yeast


[grRatio_single, grRateKO_single, grRateWT_single, hasEffect_single, delRxn_single, fluxSolution_single]...
                                   = singleGeneDeletion(ecoli_WT);                               
single_lethal_genes_ind = find(grRatio_single < 0.05 ) ; 

% Making arrays to store relevant information
result_array = zeros(length(single_lethal_genes_ind),1) ;
auxotropic_pool = cell(length(single_lethal_genes_ind),3) ; 
result_array_Gr_KO = zeros(length(single_lethal_genes_ind),1) ;
result_array_Gr_WT = zeros(length(single_lethal_genes_ind),1) ;
result_array_Ex = zeros(length(single_lethal_genes_ind),1) ;
result_array_zero_GR = zeros(length(single_lethal_genes_ind),1) ;
result_array_abundance_KO = zeros(length(single_lethal_genes_ind),1) ;
result_array_abundance_WT = zeros(length(single_lethal_genes_ind),1) ;



for i = 1:length(single_lethal_genes_ind) 
    disp(single_lethal_genes_ind(i))
    model_1 = deleteModelGenes(ecoli_WT, ecoli_WT.genes(single_lethal_genes_ind(i)));
    all_models = {model_1 , yeast_WT} ;
    [comm_models,pairedModelInfo] = create_community(all_models) ; % Creating community model for WT
    
    comm_models{1} = changeRxnBounds(comm_models{1},'EX_o2_e(u)',-0,'l'); % Set maximum oxygen uptake

    [sol_KO, result_KO] = SteadyCom(comm_models{1}) ; % Applying Steady Come to WT community model
    
    if result_KO.vBM(1) > 0 
        auxotropic_pool{i} = {comm_models{1}, model_1, result_KO} ;
        result_array(i) = 1 ;         
    end
    if result_KO.vBM(1) <= 1e-5
        result_array_zero_GR(i) = result_KO.GRmax ; 
    end
    
    GR = result_KO.vBM ;
    Ex = sum(result_KO.Ex) ; 
    %if result.GRmax == [] ; GRmax = 0; Ex = 0 ; end
    result_array_Gr_KO(i) = GR(1)/GR_WT_yeast ;
    result_array_Gr_WT(i) = GR(2)/GR_WT_ecoli ;
    result_array_Ex(i) = Ex ;
    result_array_abundance_KO(i) = result_KO.BM(1) ;
    result_array_abundance_WT(i) = result_KO.BM(2) ;
end

figure
bar(single_lethal_genes_ind, result_array_Gr_KO)
xlabel('Index of Gene deleted') ; 
ylabel('Growth Rate ratio for Mutant in Comm Model')

figure
bar(single_lethal_genes_ind, result_array_Gr_WT)
xlabel('Index of Gene deleted') ; 
ylabel('Growth Rate ratio for WT in Comm Model')

figure
bar(single_lethal_genes_ind, result_array_Ex)
xlabel('Index of Gene deleted') ; 
ylabel('Exchange Fluxes (mean)')
    
gene_ind = single_lethal_genes_ind(find(result_array)) ; 
gr_KO = result_array_Gr_KO(find(result_array)) ; 
gr_WT = result_array_Gr_WT(find(result_array)) ; 
ex_flux = result_array_Ex(find(result_array)) ;
genes = ecoli_WT.genes(gene_ind) ; 
abundance_KO = result_array_abundance_KO(find(result_array)) ;
abundance_WT = result_array_abundance_WT(find(result_array)) ;
Table = table(gene_ind,genes, gr_KO*GR_WT_ecoli, gr_WT*GR_WT_yeast,gr_KO,abundance_KO,ex_flux, 'VariableNames',...
             {'Gene Index','Gene', 'KO Growth rate' , 'WT Growth rate', 'Ratio (gr_KO/gr_WT)','KO Abundance','Sum of Ex Flux'})   

%writetable(Table,'iMM904_iML1515_result_table.csv','Delimiter',',')         
figure
bar(single_lethal_genes_ind, result_array_zero_GR)
xlabel('Index of Gene deleted') ; 
ylabel('GRmax for all combinations where GR_{KO} <= 1e-5')