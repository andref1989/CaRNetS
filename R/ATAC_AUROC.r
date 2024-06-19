#Rscript Sample_wig_file Motif_match_file Peak_file TF_bed TF_name
options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)

library(SummarizedExperiment)
library(rtracklayer)
library(igraph)
library(GenomicRanges)
source("/pbtech_mounts/khuranalab_scratch_athena/anf2034/Breast_Cancer_Network_Project/Andre_F_functions.R")
liftoverchain <- "/pbtech_mounts/khuranalab_scratch_athena/anf2034/Breast_Cancer_Network_Project/Liftover_files/hg19ToHg38.over.chain"
chain_in <- import.chain(liftoverchain)


inwig <- args[1]
motif_matches <- readRDS(args[2])
Peak_file <- readRDS(args[3])
TF_bed <- bed_to_granges(args[4])
TF_name <- args[5]
padding <-  as.numeric(args[6])
TF_bed <- liftOver(TF_bed,chain_in)

lseq <- function(from=1, to=100000, length.out=19){ exp(seq(log(from), log(to), length.out=length.out))}


#############################################################################
## Import Wig file and calculate p-values. Fixes minimum p-value at 1e-20####
#############################################################################
inwig <- import(inwig, selection=Peak_file)
inwig$pval <- 1-ppois(inwig$score,median(inwig$score))

inwig[which(inwig$pval<1e-16)]$pval <- 1e-20


#############################################################################
## Get the coordinates of peaks featuring motifs for the TF of interest######
#############################################################################
motif_ranges <- rowRanges(motif_matches)
motif_matrix <- assays(motif_matches)@listData$motifMatches

BRCA_motif_matrix <- motif_matrix[grep("BRCA_",rownames(motif_matrix)),]


BRCA_motif_ranges <- motif_ranges[grep("BRCA_",motif_ranges$name),]

index <- grep(paste0(TF_name,"_"),colnames(BRCA_motif_matrix))

if(length(index) > 1) {
    sub_mat <- BRCA_motif_matrix[,index]
    val_vec <- apply(sub_mat,1, function(x) any(x) ==TRUE)
    val_vec <- val_vec[which(val_vec ==TRUE)]
    motif_names <- rownames(sub_mat[names(val_vec),])
    str(val_vec)
    sub_ranges <- BRCA_motif_ranges[motif_names]
    print(sub_ranges)} else if (length(index) == 1){sub_mat <- BRCA_motif_matrix[,index];str(sub_mat) ;val_vec <- which(sub_mat == TRUE); motif_names <- names(sub_mat)[val_vec]; sub_ranges <- BRCA_motif_ranges[motif_names]} else {print(paste0("Cannot find motif named ",TF_name))}



###########################################################################
## AUROC for TF of interest
###########################################################################

step_counter <- lseq(1, min(inwig$pval),20)
Chip_ranges <- intersect_with_metadata(inwig,TF_bed)
ATAC_ranges <- intersect_with_metadata(inwig,sub_ranges)+padding
saveRDS(Chip_ranges,"Chip_ranges.rds")
saveRDS(ATAC_ranges,"ATAC_ranges.rds")


TP_vec <- c(1)
FP_vec <- c(1)
FN_vec <- c(0)
TN_vec <- c(0)

for(i in step_counter){
    index <- which(ATAC_ranges$pval <= i)
    if(length(index) >0){
        sub_gr <- ATAC_ranges[index]
        TP_vec <- c(TP_vec,length(intersect_with_metadata(sub_gr,Chip_ranges)))
        FP_vec <- c(FP_vec,length(setdiff_with_metadata(sub_gr,Chip_ranges)))
        FN_vec <- c(FN_vec,length(setdiff_with_metadata(Chip_ranges,sub_gr)))
        combined <- reduce(c(reduce(Chip_ranges),reduce(sub_gr)))
        TN_vec <- c(TN_vec,length(setdiff_with_metadata(inwig,combined)))}
    else {
        TP_vec <- c(TP_vec,0)
        FP_vec <- c(FP_vec,0)
        FN_vec <- c(FN_vec,length(Chip_ranges))
        TN_vec <- c(TN_vec,length(setdiff_with_metadata(inwig, Chip_ranges)))}
}

TPR <- TP_vec/(TP_vec+FN_vec)
FPR <- FP_vec/(FP_vec+TN_vec)
conting_table <- as.data.frame(cbind(TP_vec,FP_vec,FN_vec,TN_vec,TPR,FPR,TF_name),stringsAsFactors=F)

str(conting_table)

write.table(conting_table,paste0(TF_name,"_ATAC_conting_table_",padding,".txt"),row.names=FALSE,quote=F, sep='\t')
