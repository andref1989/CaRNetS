#Rscript PIQ_output_dir Fixed_output_dir promoter_file enhancer_file TF_list condensed_TFs_outdir file_extension
## Always use absolute paths

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)

PIQ_outdir <- args[1]
outdir <- args[2]
str(outdir)
outdir <- system("echo $TMPDIR",intern=TRUE)
str(outdir)

TF_list <- readLines(args[3])
condensed_outdir <- args[4]
extension <- args[5]


PIQ_outdir <- paste0(PIQ_outdir,"/")
outdir <- paste0(outdir,"/")
condensed_outdir <-paste0(condensed_outdir,"/")

ID_check <- grep("TCGA-",PIQ_outdir,value=TRUE)
if(length(ID_check) ==1){
    sample_ID <-  paste0("TCGA-",unlist(strsplit(PIQ_outdir,"TCGA-"))[2])
    sample_ID <- gsub("ATAC","ATAC_PIQ",sample_ID)} else if(length(ID_check) <1){

        int <- unlist(strsplit(PIQ_outdir,"/"))
        int <- int[length(int)]
        sample_ID <- unlist(strsplit(int,"_"))[1]
         print(sample_ID)}
 print(sample_ID)


fix_outputs <- function(indir,outdir,ext){

    system(paste0("mkdir ", outdir))

    system(paste0("cd ",indir))


    file_list <- list.files(indir)
    true_file_list <- grep(paste0("-calls.all",ext), file_list, value = TRUE)
    RC_files <- sort(grep(".RC-",true_file_list, value=TRUE))


    normal_files <- gsub(".RC-","-",RC_files)


    for (i in 1:length(normal_files)){
        TF_name = unlist(strsplit(normal_files[i],"LINE[0-9]+"))[2]
        TF_name2 = unlist(strsplit(TF_name,"-"))[1]
        TF_name = TF_name2
        print(paste0("processing ",TF_name))
        if (file.exists(paste0(outdir,TF_name,"_sorted.bed")) ==TRUE){
            print(paste0("There's already a sorted ", ext," for ",TF_name," in this directory"))}
        else{
            system(paste0("tail -n +2 ",indir, RC_files[i], " > ",indir, TF_name,"RC_int"))
            system(paste0("tail -n +2 ",indir, normal_files[i], " > ",indir, TF_name,"normal_int"))
            system(paste0("cat ", indir,TF_name,"normal_int ", indir,TF_name, "RC_int > ",indir,TF_name,"_clean.bed"))
            system(paste0("sort -V -k 1,1 -k 2,2 ", indir,TF_name,"_clean.bed > ", outdir,TF_name,"_sorted.bed"))
            system(paste0("rm ",indir,"*normal_int"))
            system(paste0("rm ",indir ,"*RC_int"))
            system(paste0("rm ",indir,"*_clean.bed"))
            quality_filter(paste0(outdir,TF_name,"_sorted.bed"),paste0(outdir,TF_name,"_qfilter.bed"))

            print(paste0("Finished processing ",TF_name))                   }}}


                                        #Cutoff set to default 700
quality_filter <- function(infile, outfile, cutoff=700){
    system(paste("awk 'FNR=1 && ($5 >= ",cutoff,")'",infile,">",outfile,sep=" "))
}

condense_motifs_to_TFs <- function(TF_list, indir,outdir){
system(paste0("mkdir ",outdir))

    suffix <- c("DN","D","I","IN")
    file_list <- list.files(indir, "_qfilter.bed")
	str(file_list)
                           
    for(i in TF_list){
        search_terms <- paste0(paste0(i,suffix),collapse="[0-9]?_|")
        in_list <- grep(search_terms,file_list,value=T)

        print(paste0("Searching for ",i,": Found ", in_list))
        if(length(in_list) == 0){
            print(paste0("No files for ",i," in this directory"))}
        else if (length(in_list) >= 1){

            system(paste0("cat ", paste0(indir,in_list,collapse=" ")," > ", paste0(outdir,i,"_condensed.bed")))}
    }}




fix_outputs(PIQ_outdir, outdir, extension)

condense_motifs_to_TFs(TF_list, outdir, condensed_outdir)
sample_ID <- gsub("/","",sample_ID)
system(paste0("tar -czvf ",PIQ_outdir, sample_ID,".gz ",PIQ_outdir,"/*.bed"))
system(paste0("rm ",PIQ_outdir,"/*-calls.all.bed"))

