#Rscript Indir Patient_specific_links Outdir

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)

indir <- paste0(args[1],"/")
p_g_links <- args[2]
outdir <- paste0(args[3],"/")
sample_id <- args[4]


reg_element_int <- function(PIQ_bed, reg_element_bed){
    #system("spack load bedtools2@2.27.0")
    check <- readLines(reg_element_bed, n=2)
    search_check <- grep("track", check)
    if(length(search_check) >0){
        new_bed <- gsub(".bed", "_fixed.bed",reg_element_bed);
        system(paste0("tail -n +2 ",PIQ_bed, " > ",new_bed))
    #    system(paste0("mv ",new_bed," ",bed))
    } else{ print("Bed file was correctly formatted")}
    system(paste0("intersectBed -a ",reg_element_bed," -b ",PIQ_bed," -wo > ", gsub("_condensed.bed","_reg_element_int.bed",PIQ_bed)))
 #   system(paste0("head ",gsub("_qfilter.bed","_reg_element_int.bed",PIQ_bed)))
}


edgelist_from_condensed <- function(indir,outfile){
    file_list <- list.files(indir,"_reg_element_int.bed")
    if(file.exists(outfile) == TRUE){
        print("The listed file already exists")}
    else{print("Creating new file")
         system(paste0("touch ", outfile))
        for(i in file_list){
        TF_name <- gsub("_reg_element_int.bed","",i)
        infile <- paste0(indir,i)
        system(paste0("awk '{print \"",TF_name,"\t\" $5 \"\t\" $4 }' ", infile," |sort > ",indir,TF_name,"_final.txt"))}
    file_list <- list.files(indir,"_final.txt")
    system(paste0("cat ",paste0(indir,file_list,collapse=" "),"|sort|uniq > ",outfile))
    system(paste0("rm ", indir,"*_final.txt"))

     }
    
            
}



file_list <- list.files(indir,"_condensed.bed")


for(i in file_list){
    index <- grep(i ,file_list)
    infile <- paste0(indir, i)

    reg_element_int(infile,p_g_links)
    print(paste0("Finished processing ",index ," of ",length(file_list)))
    
}
outfile <- paste0(outdir,sample_id,"_Final_edgelist.txt")                     
edgelist_from_condensed(indir,outfile)

system("touch Final.txt")
