# Task a.
library(dplyr) 
MakeStrangeStringList <- function(string) {
  if(!grepl("^[A-Za-z]+$", string, perl = T)) {
    cat("The input have to include only letters")
  } else {
    spliedstring <- strsplit(string, split = "") %>% 
                     unlist %>% 
                     tolower
    makeUpper <- seq(1 ,length(spliedstring), 2)
    secstring <- spliedstring
    spliedstring[makeUpper + 1] <- toupper(spliedstring[makeUpper + 1])
    secstring[makeUpper] <- toupper(secstring[makeUpper])
    return(
      list(paste(spliedstring, collapse = ""), 
           paste(secstring, collapse = "")
      )
    )
  }
}

string <- "abcdef"
MakeStrangeStringList(string) 

# Task b.
string <- "abcdef"
MakeStrangeStringList(string) 

CountLetters <- function(string) {
  sepstring <- strsplit(string, split = "") %>% 
                unlist %>% 
                as.factor
  numberofrepeatedcharacters <- length(summary(sepstring)[summary(sepstring) > 1])
  return(numberofrepeatedcharacters)
}

string <- "RhabarbArka"
CountLetters(string)


# ---TASK 1---------Preprocessing-----------------------
# This task is only a description of workflow. I haven't install those tools on my computer.
# step 1. 
# extract reads using the samtools
#   samtools required corted and indexed file:
#   samtools sort miniMNM00065.bam -o miniMNM00065.sort.bam
#   samtools index miniMNM00065.sort.bam

# and then extract reads: 
#   samtools view miniMNM00065.sort.bam chr1 -b > out.bam

# I assume that it's paired-end sequencing (no info in the instruction)
#   samtools sort -n out.bam -o out.sort.name.bam
#   samtools fastq -1 MNM00065_1.fq.gz -2 MNM00065_2.fq.gz -s /dev/null out.sort.name.bam
  
# step 2. 
# index reference chromosome using bwa and mapping using bwa mem
# for a bigger files I think is better to use mapper which return files in a .bam format because of storage limitations 
#   bwa index chr1_hg38.fasta
#   bwa mem chr1_hg38.fasta 
#         '<zcat MNM00065_1.fq.gz>' 
#         '<zcat MNM00065_2.fq.gz>' > MNM00065_aligned_pe_hg38.sam
# BWA return files in sam format. I use a samtools to convert sam to bam file.
#   samtools view -S -b MNM00065_aligned_pe_hg38.sam > MNM00065_aligned_pe_hg38.bam
