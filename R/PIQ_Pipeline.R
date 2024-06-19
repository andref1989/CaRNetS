                                        #Rscript MyScript DHS_Bamfile Motifs_file PIQ_Output_Dir Final_Output_Dir motif_symbol tmp_folder

library(stringr)
options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)


DHS_bam = args[1]
Motifs_file = args[2]
PIQ_Output_Dir = args[3]
Final_Output_Dir = args[4]
motif_symbol = args[5]
tmp = args[6]
tmpdir = system("echo $TMPDIR", intern=TRUE)
tmpdir <- paste0(tmpdir,"/")

###########
                                        #Get ID info
############
DHS_ID= grep("TCGA", DHS_bam, value = TRUE)

                                        #int <-grep("TCGA",unlist(str_split(grep("TCGA[A-Z]|[0-9]",unlist(str_split(DHS_bam, "_")), value = TRUE)[1],"[^a-zA-Z0-9']+")), value = TRUE)
if(length(DHS_ID) ==1){
int <- unlist(strsplit(DHS_ID,"TCGA-"))[2]
int <- unlist(strsplit(paste0("TCGA-",int),"\\."))[1]
DHS_ID = gsub("/","",int)
print(int)
} else if (length(DHS_ID) <1 ){
    int <- unlist(strsplit(DHS_bam,"/"))
    int <- int[length(int)]
    int <- unlist(strsplit(int,"\\."))[1]
    DHS_ID = int
    print(DHS_ID)
}

PIQ_dir <- paste0(PIQ_Output_Dir,"/",tmp,"/")
final_tmp <- tmpdir
int_dir <- paste0(tmpdir,DHS_ID,"_",tmp,"/")
print(int_dir)
system(paste0("mkdir ", PIQ_dir))
system(paste0("mkdir ", int_dir))
                                        #PIQ_funcs <- "/pbtech_mounts/khuranalab_scratch_athena/anf2034/Breast_Cancer_Network_Project/Flow/modules/PIQ/"
PIQ_funcs <- "/pbtech_mounts/homes024/anf2034/PIQ_module/"


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

###############################################################
                                        #import motif file####################
pwmtable = importJaspar(Motifs_file)

######Prep Bam File and convert to R Data format #######################################

PIQ_calls = paste0(Final_Output_Dir,"/",DHS_ID,"/")
print(PIQ_calls)

if(file.exists(paste0(PIQ_calls,DHS_ID,".RData")) == FALSE){
    system(paste0("mkdir ", Final_Output_Dir,"/",DHS_ID))
    check <- grep(".RData", DHS_bam)
    check_val <- length(check)
    if(check_val ==  0){
        print("Working 1")
        ##system(paste0("Rscript pairedbam2rdata.r common.r ",paste0(PIQ_calls,"/",DHS_ID,".RData")," ",DHS_bam))
system(paste0("Rscript bam2rdata.r common.r ",paste0(PIQ_calls,"/",DHS_ID,".RData")," ",DHS_bam))} else{ print("Working 2")
                                                                                                                 system(paste0("cp ",DHS_bam," " ,Final_Output_Dir,"/",DHS_ID))}} else { print("RData file for this TF exists in this directory")}
system(paste0("touch Finished.RData"))

system(paste0("mkdir ", final_tmp))


########################################################################################
                                        #Check motif symbol. if all, iterate through pwm file. Else, run individual motif or skip individual PWM identification.

if (motif_symbol == "all") {
    for (i in 1:length(pwmtable)){
        motif_id = i
        system(paste0("Rscript pwmmatch.exact_mod.r common.r ", paste0(Motifs_file," "), motif_id, paste0("$i ", PIQ_dir)))
##        system(paste0("Rscript pertf.r common.r ", PIQ_dir, " ",final_tmp," ", PIQ_calls," ",paste0(PIQ_calls,DHS_ID,".RData")," ",motif_id))
        system(paste0("Rscript pertf.r common.r ", PIQ_dir, " ",final_tmp," ", int_dir," ",paste0(PIQ_calls,DHS_ID,".RData")," ",motif_id))
        motif_name = unlist(strsplit(names(pwmtable[i]),"LINE[0-9]"))[2]
        unlist(strsplit(motif_name,"_"))[2]
        print(motif_name)

    }} else if (motif_symbol == "skip" ){

        for (i in 1:length(pwmtable)){
            motif_id = i
            system(paste0("Rscript pertf.r common.r ", PIQ_dir, " ",final_tmp," ", int_dir," ",paste0(PIQ_calls,DHS_ID,".RData")," ",motif_id))
##            system(paste0("Rscript pertf.r common.r ", PIQ_dir, " ",final_tmp," ", PIQ_calls," ",paste0(PIQ_calls,DHS_ID,".RData")," ",motif_id))
            motif_name = unlist(strsplit(names(pwmtable[i]),"LINE[0-9]"))[2]
            print(motif_name)
		print(motif_id)

        }} else if(motif_symbol == "none"){
                                        #TF_list <- readLines("/pbtech_mounts/khuranalab_scratch_athena/anf2034/Breast_Cancer_Network_Project/PIQ_new/Network_build_files/new_deduped_TF_list.txt")

            for (i in TF_list){

                motif_name = i
                print(motif_name)
            }} else {
                motifs_file_table = readLines(Motifs_file)
                pattern_start = grep(paste0("_",motif_symbol,"_"), motifs_file_table)
                motif_id = round((pattern_start/5)+1)
                print(motif_id)
                for(i in 1:length(motif_id)){
                        motif_sub <- motif_id[i]
print("Finding Motifs")
print(paste0("Working motif:", motif_sub))
in_cmd <- paste0("Rscript pertf.r common.r ", PIQ_dir, " ",final_tmp," ", int_dir," ",paste0(PIQ_calls,DHS_ID,".RData")," ",motif_sub)
print(in_cmd)
                                       system(paste0("Rscript pwmmatch.exact_mod.r common.r ", paste0(Motifs_file," "), motif_sub, paste0("$i ", PIQ_dir)))
##                system(paste0("Rscript pertf.r common.r ", PIQ_dir, " ",final_tmp," ", PIQ_calls," ",paste0(PIQ_calls,DHS_ID,".RData")," ",motif_id))
                
print("Footprinting now")
		
		system(in_cmd)



                motif_name = unlist(strsplit(names(pwmtable[motif_sub]),"LINE[0-9]"))[2]


                print(motif_name)
		    print(motif_id)


                    }}
print(list.files(int_dir))
system(paste0("rsync -tr ",int_dir," ",PIQ_calls))

system(paste0("rm -r ", paste0(PIQ_calls,"*.pdf")))
system(paste0("rm -r ", paste0(PIQ_calls,"*.csv")))


system("touch Finished_final.bed")

