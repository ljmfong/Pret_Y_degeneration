#### Plotting the M:F Read Depth for genes in P. reticulata ####

rm(list=ls())
ls() 

library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)

setwd("\\\\files.zoology.ubc.ca/ljmfong/flex/new_Ret_SNPs")

####### 1. Read in the file: #######
####### Look at snps first #########

#Remove genes that are putative duplicates a.k.a those that you are interested in through HI snps

gene_list_autosome <- c("ENSPREG00000000149", "ENSPREG00000000425", "ENSPREG00000000486", "ENSPREG00000020299", 
               "ENSPREG00000021037", "ENSPREG00000014816", "ENSPREG00000005754", "ENSPREG00000013757", "ENSPREG00000003677",
               "ENSPREG00000015247", "ENSPREG00000006517", "ENSPREG00000011757", "ENSPREG00000013879", "ENSPREG00000004959",
               "ENSPREG00000021017", "ENSPREG00000001126", "ENSPREG00000007572", "ENSPREG00000003854", "ENSPREG00000019543",
               "ENSPREG00000003921", "ENSPREG00000012441", "ENSPREG00000012471", "ENSPREG00000004135", "ENSPREG00000017067")
gene_list_sexchromo <- c("ENSPREG00000019446", "ENSPREG00000006420", "ENSPREG00000018782", "ENSPREG00000020252", "ENSPREG00000019378",
                         "ENSPREG00000019435", "ENSPREG00000019478", "ENSPREG00000019540", "ENSPREG00000019673",
                         "ENSPREG00000019774", "ENSPREG00000019794", "ENSPREG00000019838", "ENSPREG00000019848", "ENSPREG00000019881",
                         "ENSPREG00000019905", "ENSPREG00000019944", "ENSPREG00000019970", "ENSPREG00000019988",
                         "ENSPREG00000020032", "ENSPREG00000020049", "ENSPREG00000020126", "ENSPREG00000020131", "ENSPREG00000020223",
                         "ENSPREG00000020570", "ENSPREG00000001200", "ENSPREG00000001330", "ENSPREG00000001446", "ENSPREG00000001498", "ENSPREG00000001901", "ENSPREG00000005501",
                         "ENSPREG00000005616", "ENSPREG00000005759", "ENSPREG00000005809", "ENSPREG00000006331", "ENSPREG00000006362", "ENSPREG00000006387", "ENSPREG00000006501")

combo_gene_list <- do.call(c, list(gene_list_autosome, gene_list_sexchromo))


# Check autosomal snps and their read depth:
gene_MF_depth <- read.table(file = "check_for_dupes/MF_snps_w_genes.txt", header = TRUE)

#check male limited SNPs and if they have significantly different M_F_ratio compared to the autosomal average
LG12_male_limited <- read.table(file = "check_for_dupes/LG12_male_limited_genes.txt", header = T)
autosomal_male_snps <- read.table(file = "Autosomal_male_limited_SNPs.txt", header = T)

################################

gene_list <- c(autosomal_male_snps$gene_id, LG12_male_limited$ensembl_ID)

autosome_genes <- autosome_genes %>% 
  filter(!(geneid %in% gene_list))
mean(autosome_genes$M_F_ratio)
# 1.081317

add_chrom_column <- function(dfA, list) {
  dfA$gene_detail <- ifelse(dfA$geneid %in% list, "pute_dupe", "autosomal")
  return(dfA)
}

org_gene_MFDepth <- add_chrom_column(gene_MF_depth, gene_list)
genes_to_test <- org_gene_MFDepth %>%
  filter(geneid %in% gene_list)

autosome_genes_values <- autosome_genes$M_F_ratio
pute_dupe_values <- genes_to_test$M_F_ratio

t.test(autosome_genes_values, pute_dupe_values, alternative = "two.sided")

#Output:

data:  autosome_genes_values and pute_dupe_values
t = -2.8263, df = 64.214, p-value = 0.00627
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  -0.08939293 -0.01535638
sample estimates:
  mean of x mean of y 
1.081317  1.133691

# Any gene that has a M:F Read Depth > 1.1333691 is likely a partial dupe

write.table(genes_to_test, "check_for_dupes/organize_the_pute_dupes.txt", quote = F, row.names = F, sep = "\t")








