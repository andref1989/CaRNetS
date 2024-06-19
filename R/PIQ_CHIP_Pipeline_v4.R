#Rscript MyScript CHIP_Bamfile DHS_Bamfile Motifs_file Promoter_file Enhancer_file OutputDir PIQDir motif_symbol tmp_folder

library(stringr)
options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)

CHIP_bam = args[1] 
DHS_bam = args[2]
Motifs_file = args[3]
promoter_bed = args[4]
enhancer_bed = args[5]
Output_dir = args[6]
PIQ_dir = args[7]
motif_symbol = args[8]
tmp = args[9]


Output_dir <- paste0(Output_dir,"/")
PIQ_dir <- paste0(PIQ_dir,"/",tmp,"/")
system(paste0("mkdir ", PIQ_dir))
PIQ_funcs <- "/pbtech_mounts/khuranalab_scratch_athena/anf2034/Breast_Cancer_Network_Project/Flow/modules/PIQ/"


system(paste0("cp ", PIQ_funcs,"*.r ."))
##Define important functions

##Import PWMs
importJaspar <- function(file=myloc) {
      vec <- readLines(file)
        vec <- gsub("\t"," ",vec)
        vec <- gsub("\\[|\\]", "", vec)
        start <- grep(">", vec); end <- grep(">", vec) - 1
        pos <- data.frame(start=start, end=c(end[-1], length(vec)))
        pwm <- sapply(seq(along=pos[,1]), function(x) vec[pos[x,1]:pos[x,2]])
        pwm <- sapply(seq(along=pwm), function(x) strsplit(pwm[[x]], " {1,}"))
        pwm <- lapply(seq(along=start), function(x) matrix(as.numeric(t(as.data.frame(pwm[(pos[x,1]+1):pos[x,2]]))[,-1]), nrow=4, dimnames=list(c("A", "C", "G", "T"), NULL)))
        names(pwm) <- gsub(">", "", vec[start])
        return(pwm)
  }


##import PIQ calls and push out list of genes only
build_gene_list <- function(motif_symbol,dir,output_name) {
    TF_calls_bed = grep(motif_symbol, list.files(dir), ignore.case = TRUE,value=TRUE)
    PIQ_calls_bed = grep("-calls.all.bed", TF_calls_bed, ignore.case=TRUE, value=TRUE)
    if(length(PIQ_calls_bed) == 0) { print("No PIQ calls for this TF in specified directory") } else if (length(PIQ_calls_bed) == 1){ system(paste0("tail -n +2 ",dir, PIQ_calls_bed, " > " , paste0(dir,motif_symbol,"_clean.bed")))} else { for (i in 1:length(PIQ_calls_bed)) { system(paste0("tail -n +2 ",dir,PIQ_calls_bed[i]," >> ", paste0(dir,motif_symbol,"_clean.bed")))}}
    system(paste0("sort -V -k 1,1 -k 2,2 ",dir, motif_symbol,"_clean.bed > ",dir, motif_symbol,"_all_motif_calls_final.bed"))
     
    
}


Concatenate_final <- function(TF_list, input_dir, output_dir){
    TF_list2 = readLines(TF_list)
    for (i in TF_list2){
        system(paste0("cat ",paste0(input_dir,i,"_enhancer_genes_final.bed "), paste0(input_dir,i,"_promoter_genes_fina
l.bed > "), paste0(output_dir,i,"_combined_regelement.bed")))}}


reduce_TF_list_func <- function (TF_list, input_dir, output_dir){
    TF_list2 = readLines(TF_list)
    for (i in TF_list2){
        system(paste0("cat ",input_dir,i,"_combined_regelement.bed |sort|uniq > ", output_dir,i, "_reduced_combo_regelement.bed"))
        system(paste0("awk -v OFS='\t' '{print $2, $3}' ",output_dir,i,"_reduced_combo_regelement.bed |sort|uniq > ",output_dir,i,"_no_regelementinfo.bed"))

    }}


###############################################################
#import motif file####################
pwmtable = importJaspar(Motifs_file)

######Prep Bam File and convert to R Data format #######################################
DHS_ID= grep("GSM", DHS_bam, value = TRUE)
int <-grep("GSM",unlist(str_split(grep("GSM[0-9]",unlist(str_split(DHS_bam, "_")), value = TRUE)[1],"[^a-zA-Z0-9']+")), value = TRUE)
DHS_ID = int
PIQ_calls = paste0(Output_dir,DHS_ID,"/")

if(file.exists(paste0(Output_dir,DHS_ID,"/",DHS_ID,".RData")) == FALSE){
    system(paste0("Rscript bam2rdata.r common.r ",paste0(Output_dir,DHS_ID,"/",DHS_ID,".RData")," ",DHS_bam))} else { print("RData file for this TF exists in this directory")}

system(paste0("mkdir ", tmp))
######format and import promoter file and enhancer file#############
system(paste0("mkdir ", Output_dir,"/",DHS_ID))
if(file.exists(paste0(PIQ_calls,"promoter.bed")) == FALSE) {
system(paste0("awk -v OFS='\t' '{print $1, $2, $3, $4}' ", promoter_bed, " > ",PIQ_calls , "promoter.bed"))} else {file.info(paste0(PIQ_calls,"promoter.bed"))}
if(file.exists(paste0(PIQ_calls,"enhancer.bed")) == FALSE) {
system(paste0("awk -v OFS='\t' '{print $1, $2, $3, $4}' ", enhancer_bed , " > ",PIQ_calls , "enhancer.bed"))} else {file.info(paste0(PIQ_calls,"enhancer.bed"))}



########################################################################################
#Check motif symbol. if all, iterate through pwm file. Else, run individual motif or skip individual PWM identification.
if (motif_symbol == "all") {
    for (i in 1:length(pwmtable)){
        motif_id = i
        system(paste0("Rscript pwmmatch.exact.r common.r ", paste0(Motifs_file," "), motif_id, paste0("$i ", PIQ_dir)))
        system(paste0("Rscript pertf.r common.r ", PIQ_dir, " ",tmp," ", PIQ_calls," ",paste0(PIQ_calls,DHS_ID,".RData")," ",motif_id))
        motif_name = unlist(strsplit(names(pwmtable[i]),";"))[2]
        print(motif_name)
#        build_gene_list(motif_name,PIQ_calls,paste0("/pbtech_mounts/khuranalab_scratch_athena/anf2034/Breast_Cancer_Network_Project/PIQ_new/Gene_Lists/",DHS_ID))
        }} else if (motif_symbol == "skip" ){

           for (i in 1:length(pwmtable)){
            motif_id = i
        system(paste0("Rscript pertf.r common.r ", PIQ_dir, " ",tmp," ", paste0(Output_dir,DHS_ID,"/ "), paste0(Output_dir,DHS_ID,"/",DHS_ID,".RData")," ",motif_id))
            motif_name = unlist(strsplit(names(pwmtable[i]),";"))[2]
            print(motif_name)
            build_gene_list(motif_name,PIQ_calls,paste0("/pbtech_mounts/khuranalab_scratch_athena/anf2034/Breast_Cancer_Network_Project/PIQ_new/Gene_Lists/",DHS_ID))

        }} else if(motif_symbol == "none"){
            TF_list <- readLines("/pbtech_mounts/khuranalab_scratch_athena/anf2034/Breast_Cancer_Network_Project/PIQ_new/Network_build_files/new_deduped_TF_list.txt")

            for (i in TF_list){
            
            motif_name = i
            print(motif_name)
            build_gene_list(motif_name,"/pbtech_mounts/khuranalab_scratch_athena/anf2034/Breast_Cancer_Network_Project/PIQ_new/ESR_run/GSM1255280/","/pbtech_mounts/khuranalab_scratch_athena/anf2034/Breast_Cancer_Network_Project/PIQ_new/Gene_Lists/GSM1255280_test/")}} else {
            motifs_file_table = readLines(Motifs_file)
            pattern_start = grep(motif_symbol, motifs_file_table)[1]
            motif_id = round((pattern_start/5)+1)
    system(paste0("Rscript pwmmatch.exact.r common.r ", paste0(Motifs_file," "), motif_id, paste0("$i ", PIQ_dir)))
        system(paste0("Rscript pertf.r common.r ", PIQ_dir, " ",tmp," ", paste0(Output_dir,DHS_ID,"/ "), paste0(Output_dir,DHS_ID,"/",DHS_ID,".RData")," ",motif_id))


          motif_name = unlist(strsplit(names(pwmtable[i]),";"))[2]
            print(motif_name)
            build_gene_list(motif_name,PIQ_calls,paste0("/pbtech_mounts/khuranalab_scratch_athena/anf2034/Breast_Cancer_Network_Project/PIQ_new/Gene_Lists/",DHS_ID))

        }



