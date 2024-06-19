#Rscript Sample_wig_file Motif_match_file Peak_file TF_list Peak_Gene_links
## Always use absolute paths
start <- Sys.time()
options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)

library(SummarizedExperiment)
library(rtracklayer)
library(igraph)
source("/pbtech_mounts/khuranalab_scratch_athena/anf2034/Breast_Cancer_Network_Project/Andre_F_functions.R")

inwig <- args[1]
motif_matches <- readRDS(args[2])
Peak_file <- readRDS(args[3])
TF_list <- sort(readLines(args[4]))
Peak_Gene_links <- readRDS(args[5])
outdir <- args[6]

if(file.exists(outdir) == FALSE){
    system(paste0("mkdir ",outdir))} else{ print("Output directory already exists")}
filename <- inwig
inwig <- import(inwig,selection=Peak_file)
inwig$pval <- 1-ppois(inwig$score,median(inwig$score))
inwig_sig <- inwig[which(inwig$pval < 1e-10)]

motif_ranges <- rowRanges(motif_matches)
motif_matrix <- assays(motif_matches)@listData$motifMatches

BRCA_motif_matrix <- motif_matrix[grep("BRCA_",rownames(motif_matrix)),]


BRCA_motif_ranges <- motif_ranges[grep("BRCA_",motif_ranges$name),]

edgelist <- data.frame(stringsAsFactors=F)
for(i in TF_list){
    print(paste0("Working ", i))
    sub_mat <- BRCA_motif_matrix[,grep(paste0(i,"_"),colnames(BRCA_motif_matrix))]
    
    if(is.null(dim(sub_mat)) == FALSE){
        print(paste0("Multiple motifs for ",i))

        val_vec <- apply(sub_mat,1, function(x) any(x) ==TRUE)
        
        val_vec <- val_vec[which(val_vec ==TRUE)]
        motif_names <- rownames(sub_mat)[val_vec]}
 
   
    else { val_vec <- which(sub_mat == TRUE)
           motif_names <- names(sub_mat)[val_vec]
       }



    sub_ranges <- BRCA_motif_ranges[motif_names,]
    sub_ranges2 <- intersect_with_metadata(sub_ranges,inwig_sig)

    range_int <- intersect_with_metadata(Peak_Gene_links,sub_ranges2)

    reg_targets <- unique(sort(range_int$Target))
    e_sub <- data.frame(cbind(i,reg_targets))
    if(dim(e_sub)[2] == 2){
        edgelist <- rbind(edgelist,e_sub)}
    else if (dim(e_sub)[2] == 1){
        print(paste0("No regulatory targets for ",i))}
    print(paste0("Finished Processing ", i))
}


outfile <- unlist(strsplit(filename,"Bigwigs/"))[2]
outfile <- gsub(".bw","",outfile)

sample_graph <- graph_from_edgelist(as.matrix(edgelist))
#write.table(edgelist,paste0(outfile,"_edgelist.txt"),quote=F,row.names=F,col.names=F,sep='\t')

saveRDS(sample_graph,paste0(outdir,outfile,"_graph.rds"))
end <- Sys.time()
elapsed <- end-start
print(elapsed)
    
