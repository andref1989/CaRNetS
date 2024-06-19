##MyScript Achilles_data genes_of_interest DGI_df TF_score_df CGC_list TF_list subset corr_cutoff output_dir Achilles_metadata Expression_df


source("~/Andre_F_functions.R")
require(dplyr)


options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)


Achilles_data <- readRDS(args[1])
genes_of_interest <- readLines(args[2])
DGI_df <- readRDS(args[3])
TF_score_df <- readRDS(args[4])
CGC_list <- readLines(args[5])
TF_list <- readLines(args[6])
subset <- args[7]
corr_cutoff <- as.numeric(args[8])
outdir <- paste0(args[9],"/")
Achilles_metadata <- readRDS(args[10])
Expression_df <- readRDS(args[11])




if(is.na(file.info(outdir)$isdir) == TRUE){ system(paste0("mkdir ", outdir))} else{ print(paste0("Outputting file to ", outdir))}

#######################################
## Setting fixed parameters for NTP####
#######################################
TF_score_column <- "Disease"
template_length <- 150
most_variable_gene_cutoff <- 0.75
random_sampling <- 1000
verbose=TRUE
num_predictions_filter <- NULL
cor_method <- "spearman"
output_length <- 3000
##NTP_file <- "/pbtech_mounts/homes024/anf2034/NTP_results.rds"
##NTP_file <- "/pbtech_mounts/homes024/anf2034/NTP_results_6_16_final.rds"
NTP_file <- "/pbtech_mounts/homes024/anf2034/NTP_results_7_2_final.rds"
########################################


if(file.exists("/pbtech_mounts/homes024/anf2034/NTP_results.rds")==FALSE){
NTP_all <- nearest_template_prediction_v2(Expression_df,Expression_df[setdiff(colnames(Expression_df), TF_score_df$Sample)], TF_score_df,cluster_column=TF_score_column,is_ranked=FALSE,method="simple",template_length=template_length,quantile_cutoff=most_variable_gene_cutoff,random_sample=random_sampling,verbose=verbose)[[2]]
NTP_results <- NTP_all[1]
saveRDS(NTP_all,"NTP_results.rds")
} else{
print("Importing NTP results")
##NTP_results <- readRDS("NTP_results.rds")[1]
}

##str(NTP_results)


##NTP_results_final <- parse_NTP(NTP_results,Achilles_metadata,0.1)[[2]]
##NTP_results_final <- filter(NTP_results_final, Signif==TRUE)
##str(NTP_results_final)
NTP_results_final <- readRDS(NTP_file)
NTP_results_final <- filter(NTP_results_final, Signif==TRUE)
table(NTP_results_final$Prediction)

if(is.null(num_predictions_filter)==TRUE){ NTP_results_final <- NTP_results_final} else { NTP_results_final <- filter(NTP_results_final, Num_Predictions >= num_predictions_filter)}


coessential_genes <- get_coessential_genes_v2(Achilles_data, NTP_results_final, subset=subset, genes_of_interest, corr_cutoff=corr_cutoff,cor_method=cor_method)
druggable_genes <- get_druggable_genes(coessential_genes,DGI_df,Achilles_data)
druggable_df <- get_final_candidate_scores_v2(druggable_genes[[1]], TF_score_df, Achilles_data, CGC_list, TF_list,subset=subset,flag=coessential_genes[[5]],cell_line_mapping_df=NTP_results_final)

if(output_length=="All"){
write.table(druggable_df,paste0(outdir,subset,"_druggable_drivers_and_coessential_gene_candidates.txt"), quote=F, row.names=F, sep='\t')} else{
write.table(druggable_df[1:max(output_length,nrow(druggable_df)),],paste0(outdir,subset,"_druggable_drivers_and_coessential_gene_candidates.txt"), quote=F, row.names=F, sep='\t')}

print("All Done")

