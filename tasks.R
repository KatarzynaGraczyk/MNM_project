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

# -------- TASK2 ------- Structural variants hands-on --------------
library(vcfR)
library(ggplot2)
library(stringr)

# read vcf file 
vcf <- read.vcfR( "tumor_vs_normal.manta.somatic.vcf.gz", verbose = FALSE )
vcfFix <- getFIX(vcf, getINFO = T) %>% as.data.frame()

# subtask 1
numer_of_breakends <- sum(grepl("SVTYPE=BND", vcfFix$INFO))/2

# subtask 2
sel_deletions <- vcfFix[grepl("VTYPE=DEL", vcfFix$INFO), ]
sel_deletions$del_len <- sel_deletions$REF %>% nchar()
sel_deletions$CHROM <- factor(sel_deletions$CHROM, levels = unique(sel_deletions$CHROM)) 

ggplot(sel_deletions, aes(x = CHROM, y = del_len)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.3), size = 0.3, color = "darkgreen") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Boxplot of the deletion length per chromosome", 
       x = "chromosome", y = "deletion length")
  
# subtask 3
snp_fail <- vcfFix[vcfFix$FILTER != "PASS", ] 
cat(nrow(snp_fail), "variants failed to pass the filtering")

snp_fail_stats <- table(snp_fail[, "FILTER"]) %>% 
                    as.data.frame()

ggplot(snp_fail_stats, aes(x = "", y=Freq, fill=Var1)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + 
  theme(legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5)) + 
  labs(fill = "Causes of failing:") +
  scale_fill_manual(values = c("#d8b365", "#999999", "#5ab4ac")) + 
  labs(title="Piechart of the most frequent reasons to fall the filtering") 

# subtask 4
cipos <- vcfFix[grepl("CIPOS", vcfFix$INFO), ]
cipos$CI <- str_extract(cipos$INFO, "CIPOS=.[0-9]*") %>% 
  sub("CIPOS=", "", .) %>%
  as.numeric()
sel_variant <- cipos[cipos$CI == min(cipos$CI), ]

# subtask 5
vcfFix[vcfFix$ID == "MantaBND:28842:0:1:0:0:0:0", ]
# type of variant: breakends, reverse complement piece extending left of X:20306777 is joined after T
