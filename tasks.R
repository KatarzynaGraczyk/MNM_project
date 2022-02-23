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


# Task b.
string <- "abcdef"
MakeStrangeStringList(string) 

CountLetters <- function(string) {
  sepstring <- strsplit(string, split = "") %>% unlist %>% as.factor
  numberofrepeatedcharacters <- length(summary(sepstring)[summary(sepstring) > 1])
  return(numberofrepeatedcharacters)
}

string <- "RhabarbArka"
CountLetters(string)


# ---TASK 1---------Preprocessing-----------------------
# step 1. 
# using a samtools extract reads 
# samtools required corted and indexed file:
#samtools sort miniMNM00065.bam -o miniMNM00065.sort.bam
#samtools index miniMNM00065.sort.bam

# and then extract reads: 
#samtools view miniMNM00065.sort.bam chr1 -b > out.bam

# I assume that it's paired-end sequencing (no info in the instruction)
#samtools sort -n out.bam -o out.sort.name.bam
#samtools fastq -1 MNM00065_1.fq.gz -2 MNM00065_2.fq.gz -s /dev/null out.sort.name.bam
  
# step 2. 
# index reference chromosome using bwa and mapping using bwa mem. 
# For a bigger files I think is better to use mapper which return files in a .bam format because of storage limitations. 

#bwa index chr1_hg38.fasta
#bwa mem chr1_hg38.fasta 
#         '<zcat MNM00065_1.fq.gz>' 
#         '<zcat MNM00065_2.fq.gz>' > MNM00065_aligned_pe_hg38.sam
# BWA return files in sam format. I use a samtools to convert sam to bam file.
#samtools view -S -b MNM00065_aligned_pe_hg38.sam > MNM00065_aligned_pe_hg38.bam


# ------ TASK3 ------ Structural variants hands-on -----------------
library(vcfR)
vcf <- read.vcfR( "tumor_vs_normal.strelka.somatic.snvs.vcf.gz", verbose = FALSE )
head(vcf)

# 1) filter reads ...
snp_pass <- vcf@fix[, 7] == "PASS"
selected_variants <- vcf[snp_pass, ]
selected_variants_to_annotations <- selected_variants@fix

#2) make an anotation using biomart
library(biomaRt)
ensembl_snp=useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")

chr_name <- selected_variants_to_annotations[,1] %>% as.numeric()
allele <-  selected_variants_to_annotations[,4]
motif_start <- selected_variants_to_annotations[,2] %>% as.numeric()

annotations <- getBM(attributes=c("refsnp_source",
                   'refsnp_id',
                   'chr_name',
                   'chrom_start',
                   'chrom_end',
                   "consequence_type_tv", 
                   "clinical_significance"), 
      values = c(chr_name, allele, motif_start),
      mart = ensembl_snp)






