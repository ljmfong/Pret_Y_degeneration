##########################################
########### Amino Acid Changes ###########
##########################################


  # 1. Load libraries & Main File

library(stringr)
library(dplyr)
library(tidyverse)
library(gsubfn)
library(ggplot2)
library(RColorBrewer)

rstudioapi::writeRStudioPreference("data_viewer_max_columns", 1000L)

setwd("//files.zoology.ubc.ca/ljmfong/flex/new_Ret_SNPs")
getwd()

amino_acid_key <- read.table(file.choose(), sep = "\t", header = T)
#protein_check_for_SNPs/aa_change_key.txt

init_snp_change <- read.table(file.choose(), sep = "\t", header = T)
#protein_check_fo_SNPs/version_1_SNP_aa_change.txt

count_aa_change <- init_snp_change %>%
  inner_join(amino_acid_key, by = c("aa_org", "aa_chg"))
count_aa_change$sum <- (count_aa_change$structure) + (count_aa_change$polarity) + (count_aa_change$charge)

# Remove sex-specific SNPs:

autosome_aa_change <- subset(count_aa_change, CHROM != "LG12") 

# The sex-specific SNPs on autosomes + duplicates that need to be removed were found at: 
# LG4	20150922
# LG5	3369885
# LG7	2860383, 2860420, 2860509, 2860512, 2860711, 2860642, 2860921
# LG8 22642317
# LG22  7652342

other_sexlim_list <- c(20150922, 3369885, 2860383, 2860420, 2860509, 2860512, 2860711, 2860642, 2860921, 22642317, 7652342)

autosome_aa_change <- autosome_aa_change %>%
  filter(!POS %in% other_sexlim_list) 

sexchromo_aa_change <- subset(count_aa_change, CHROM == "LG12")

sexlim_list <- c(5652422, 5705274, 5705331, 20712005, 20931822, 20936465, 20936621, 20938276, 20942935, 21072080, 21214735, 21285375,
              21285399, 21328421, 21340862, 21344523, 21344563, 21345060, 21348667, 21359254, 21364902, 21365941, 21372888,
              21385522, 21385618, 21385891, 21387147, 21387430, 21387661, 21391427, 21391639, 21438860, 21443385, 21466508,
              21466550, 21468271, 21468532, 21486526, 21513934, 21514479, 21514515, 21515845, 21888711, 24831338, 24923220,
              24964412, 24964761, 24996712, 24998564, 25141492, 25141630, 25141935, 25141951, 25261456, 25261461, 25267942,
              25275647, 25278301, 25280384, 25280689, 25281082, 25281125, 25281231, 25281902, 25281946, 25319803, 25319966,
              25321077, 25833552, 25836614, 25840579, 25861233, 25863124)

sexchromo_aa_change <- sexchromo_aa_change %>%
  filter(!POS %in% sexlim_list) #I have 5,394 SNPs


nonsex_specific_snps <- rbind(autosome_aa_change, sexchromo_aa_change)

    ### Look at the percent change in the amino acids - overall ###

  # Non-sex specific:

count(nonsex_specific_snps, sum == 0) # 31061 show no change across all properties 
count(nonsex_specific_snps, sum == 1) # 45645 show one change across all properties
count(nonsex_specific_snps, sum == 2) # 21766 show two changes across all properties 
count(nonsex_specific_snps, sum == 3) # 2688 show change across all properties

  # Compared to your missense SNPs that are sex-limited (66):

# 0 changes = 20/66
# 1 change = 29/66
# 2 changes = 12/66 
# 3 changes = 5/66 

# Structure = 19/66 
# Polarity = 31/66 
# Charge = 19/66

#### Calculate Stats ###

snp_table <- data.frame(category = c("Non-sex-limited", "Male-limited"),
                        Silent = c(140422,83),
                        Missense = c(101160,66),
                        LOF = c(3748,4))
snp_rows <- snp_table[,-1]
rownames(snp_rows) <- snp_table[,1]

chisq.test(snp_rows)
#Output:

#Pearson's Chi-squared test

#data:  snp_rows
#X-squared = 1.5543, df = 2, p-value = 0.4597

fisher.test(snp_rows, simulate.p.value = TRUE, B = 1e5) # p-value = 0.3823






