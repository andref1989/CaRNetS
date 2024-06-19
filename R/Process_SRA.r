                                        #Rscript MyScript SRA_ID GEO_ID Sample_info(cell type) genome build
options(echo=TRUE)
args <- commandArgs(trailingOnly=TRUE)
print(args)

SRA_ID <- args[1]
GEO_ID <- args[2]
Sample_info <- args[3]
genome_build <- args[4]
outdir <- args[5]
Bowtie_index_dir <- "/athena/khuranalab/scratch/anf2034/Breast_Cancer_Network_Project/BOWTIE2_INDEXES/"


                                        #system("spack versions bowtie2")
                                        #system("spack add bowtie2")
                                        #system("spack commands")
system("source ~/.bash_profile")
system("which bowtie2")
system("which fastqc")
system("which trim_galore")

                                        #Download and pipe file to fastq

in_fastq <- paste0(outdir,GEO_ID,"_",Sample_info,"_DNase_",genome_build,".fastq")
out_fastq <- paste0(outdir,GEO_ID,"_",Sample_info,"_DNase_",genome_build,"_trimmed.fastq")
if( file.exists(in_fastq) == FALSE){
    system(paste0("fastq-dump -O ",outdir," ",SRA_ID))
    system(paste0("mv ", outdir, SRA_ID,".fastq ", in_fastq))
} else {print("File already exists")}


## Trim adapter sequences from FASTQ file
if (file.exists(out_fastq) == FALSE){
    int_dir <- paste0(outdir,"/",Sample_info,"_",GEO_ID,"/")
    system(paste0("mkdir ", int_dir))
    file_list <- list.files(int_dir,".fq$")
    str(file_list)
    if(!file.exists(paste0(int_dir,file_list[1]))){
        print("Trimming")
        system(paste0("trim_galore ", in_fastq," -o ",int_dir))} else{
                
            print(paste0("Using ",file_list))}

    file_list <- list.files(int_dir,".fq$")
    system(paste0("mv ", int_dir,file_list[1]," ",in_fastq))
           
    system(paste0("fastqc ", in_fastq))
    qc_file <- gsub(".fastq","_fastqc.zip",in_fastq)

    system(paste0("unzip ",qc_file," -d ",outdir))
    qc_out <- readLines(paste0(gsub(".zip","",qc_file),"/fastqc_data.txt"))

    line1 <- grep("Overrepresented",qc_out)
    line2 <- grep("END_MODULE", qc_out)
    diff <- line2-line1
    line2 <- min(line2[diff > 0])
    qc_sub <- qc_out[(line1+1):(line2-1)]
    writeLines(gsub("#","",qc_sub),"qc_sub")
    qc_df <- read.table("qc_sub", sep='\t',header=T)
    a_seq <- as.character(qc_df[,1])
    print(a_seq)
    if(a_seq == 'pass'){
        print("No overrepresented sequences detected")
        system(paste0("mv ", in_fastq," ",out_fastq))
        system(paste0('echo  "Processed" > ', in_fastq))
    } else if (a_seq != 'pass'){
        adapters <- paste0("-b ", qc_df[,1], collapse= ' ')
        system(paste("cutadapt", adapters,"-o",out_fastq,in_fastq,collapse=" "))}} else {print(paste0("Adapters have already been trimmed for this file, please check ", out_fastq))}

dir <- paste0(getwd(),"/")
## Align Fastq file with Bowtie to reference genome in bowtie index dir

system(paste0("cd ",Bowtie_index_dir))
out_sam <- paste0(outdir,GEO_ID,"_",Sample_info,"_",genome_build,".sam")

if(file.exists(out_sam)== FALSE){
    system(paste0("bowtie2 -x ",Bowtie_index_dir,genome_build," --very-fast-local -U ",out_fastq," -S ", out_sam))} else{ print(paste0("File with name ", out_sam," already exists"))}
if(file.exists(gsub(".sam",".bam",out_sam)) == FALSE){
    system(paste0("samtools view -bS ",out_sam," > ",gsub(".sam",".bam",out_sam)))} else{ print(paste0("File with name ", gsub(".sam",".bam",out_sam)," already exists"))}


                                        #system(paste0("gzip ",dir,out_fastq))

