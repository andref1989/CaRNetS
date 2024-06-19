#Rscript TF Indir Outdir

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)

TF <- args[1]
indir <- paste0(args[2],"/")
outdir <- paste0(args[3],"/")

library(igraph)

make_TF_target_mat <- function(TF,indir="/pbtech_mounts/homes024/anf2034/ATAC_Seq_Project/patientEdgeList/",outdir="/pbtech_mounts/homes024/anf2034/ATAC_Seq_Project/Target_matrices/"){
    require(igraph)
    source("/home/anf2034/Andre_F_functions.R")

    target_list <- list()
    print(paste0("Working on ", TF))

file_list <- list.files(indir,"_graph.rds")
#file_list <- file_list[grep("_graph",file_list)]
#file_list <- file_list[1:10]

for(i in file_list){graph <- readRDS(paste0(indir,i))
                    targets <- get_regulatory_targets(graph,TF)
                    target_list <- c(target_list,list(targets))
                    print(paste0("Finished ", grep(i,file_list)," of ",length(file_list)))

                }
saveRDS(target_list, paste0(outdir,TF,"_target_list.rds"))


names(target_list) <- gsub("_graph.rds","", file_list)

union_targets <- c();  for(i in target_list){ union_targets <- unique(c(union_targets,i))}
TF_matrix <- matrix(0,length(target_list),length(union_targets))
rownames(TF_matrix) <- names(target_list)
colnames(TF_matrix) <- union_targets
str(TF_matrix)
for(i in 1:length(target_list)){ targets <- target_list[[i]]
	 if(length(targets) >1) { TF_matrix[i, targets] <- 1 } else{ TF_matrix[i,] <- 0}}

target_summary <- colSums(TF_matrix)/nrow(TF_matrix)

print(paste0("Variance for ", TF,"= ",round(var(target_summary),digits=2)))

saveRDS(TF_matrix,paste0(outdir,TF,"_target_matrix.rds"))
saveRDS(target_summary,paste0(outdir,TF,"_target_summary.rds"))


}

make_TF_target_mat(TF, indir,outdir)
