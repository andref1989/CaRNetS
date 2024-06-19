#Rscript enhancer_bed promoter_bed CHIP_TF CHIP_file PIQ_output_dir output_dir

library(stringr)
library(pracma)


#Load Arguments
options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)

enhancer_bed <- args[1]
promoter_bed <- args[2]
CHIP_TF <- args[3]
CHIP_file <- args[4]
PIQ_output_dir <- args[5]
PIQ_output_dir <- paste0(PIQ_output_dir,"/")
output_dir <- args[6]

outfile <- paste0(output_dir,"combined_regelement.bed")
system(paste0("mkdir ",output_dir))

make_reg_element_file <- function(enhancer_bed, promoter_bed, output_dir){
    outfile <- paste0(output_dir,"combined_regelement.bed")
    system(paste0("cat ",enhancer_bed," ",promoter_bed," > ", outfile))
    system(paste0("sortBed -i ", outfile," > ",output_dir,"sorted_outfile"))
    system(paste0("mv ",output_dir,"sorted_outfile ", outfile))
    }

make_gene_connections <- function(CHIP_file, CHIP_TF,PIQ_output_dir, output_dir){
    outfile <- paste0(output_dir,"combined_regelement.bed")
    CHIP_out <- paste0(output_dir,CHIP_TF,"_CHIP.bed")
    PIQ_out <- paste0(output_dir,CHIP_TF,"_PIQ.bed")
    system(paste0("sortBed -i ",CHIP_file," > ", CHIP_out))
    system(paste0("sortBed -i ",PIQ_output_dir, CHIP_TF,"_condensed.bed > ",PIQ_out))
    system(paste0("bedtools intersect -a ", outfile, " -b ", CHIP_out, " > ", output_dir,CHIP_TF,"_chip_sites_genes.bed"))
    system(paste0("bedtools intersect -a ", outfile, " -b ", PIQ_out, " -wo > ", output_dir,CHIP_TF,"_PIQ_sites_genes.bed"))
}

plot_AUROC <- function(output_dir, CHIP_TF){
    pdf(paste0(output_dir,CHIP_TF,"_AUROC.pdf"))

   CHIP_in <- read.table(paste0(output_dir,CHIP_TF,"_chip_sites_genes.bed"), sep = "\t", stringsAsFactors=FALSE)
   PIQ_in <- read.table(paste0(output_dir,CHIP_TF,"_PIQ_sites_genes.bed"), sep = "\t", stringsAsFactors=FALSE)
   

   CHIP_in_file <- paste0(output_dir,CHIP_TF,"_chip_sites_genes.bed")

   cutoff_min <- min(PIQ_in$V9)
   cutoff_max <- max(PIQ_in$V9)
   AUC_step_counter <- seq(cutoff_min, cutoff_max, cutoff_max/250)

   system(paste0("awk '{print $4}' ",CHIP_in_file," |sort|uniq > ",output_dir, CHIP_TF, "_Chip_Genes"))
   system (paste0("awk '{print ","\"",CHIP_TF,"\t\" $1}' ",output_dir, CHIP_TF, "_Chip_Genes > ", output_dir ,CHIP_TF,"_Chip_Genes.txt"))
   CHIP_Seq_truthset = read.table(paste0(output_dir,CHIP_TF,"_Chip_Genes.txt"), sep = "\t", stringsAsFactors = FALSE)
   CHIP_genes <- CHIP_Seq_truthset$V2

  TP_VEC <- length(CHIP_genes)
  FP_VEC <- length(CHIP_genes)
  TN_VEC <- 0
  FN_VEC <- 0
  PIQ_min <-0
  PIQ_max <- 1000

   for (i in AUC_step_counter) {
    Subset_PIQ_predictions = subset(PIQ_in, PIQ_in$V9 >= i)
    Subset_genes <- unique(Subset_PIQ_predictions$V4)
    common_genes <- intersect(Subset_genes,CHIP_genes)

    TP <- length(common_genes)
    FP <- length(Subset_genes)-length(common_genes)
    FN <- length(CHIP_genes)- TP
    TN <- 20246 - (TP+FP+FN)
    min_val <- i
    max_val <- 1000

    TP_VEC <- c(TP_VEC, TP)
    FP_VEC <- c(FP_VEC, FP)
    FN_VEC <- c(FN_VEC, FN)
    TN_VEC <- c(TN_VEC, TN)
    PIQ_min <- c(PIQ_min, min_val)
    PIQ_max <- c(PIQ_max , max_val)

    print(i)}

    TPR_VEC <- TP_VEC/(TP_VEC+FN_VEC)
    FPR_VEC <- FP_VEC/(FP_VEC+TN_VEC)
       
   contingency_table <- cbind(TP_VEC,FP_VEC,TN_VEC,FN_VEC,TPR_VEC,FPR_VEC, PIQ_min, PIQ_max)

   write.table(contingency_table , paste0(output_dir,CHIP_TF,"_contingency_table.txt"), row.names=FALSE, quote=FALSE)


   TPR_VEC <- c(1,TPR_VEC,0)
   FPR_VEC <- c(1,FPR_VEC,0)


   x = seq(0,1,0.1)
   y = seq(0,1,0.1)

   AUC_PRO_ENH <- abs(trapz(FPR_VEC, TPR_VEC))
   AUC_PRO_ENH <- round(AUC_PRO_ENH, digits = 3)

   cutoff_diff <- abs(PIQ_min - 700)
   cutoff_index <- which(cutoff_diff == min(cutoff_diff))

   plot(FPR_VEC, TPR_VEC, xlim = c(0,1), ylim = c(0,1), main = paste0(CHIP_TF,":AUC =", AUC_PRO_ENH), xlab = "FPR", ylab="TPR",type = "l")
lines (x,y, lty = 2, col = "blue")

}


make_reg_element_file(enhancer_bed, promoter_bed, output_dir)
make_gene_connections(CHIP_file, CHIP_TF,PIQ_output_dir, output_dir)
plot_AUROC(output_dir, CHIP_TF)
