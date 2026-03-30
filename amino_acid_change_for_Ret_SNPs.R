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
  filter(!POS %in% other_sexlim_list) #I have 95,693 SNPs

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


nonsex_specific_snps <- rbind(autosome_aa_change, sexchromo_aa_change) # I have 101,160 SNPs

    ### Look at the percent change in the amino acids - overall ###

  # 0. Non-sex specific:

count(nonsex_specific_snps, sum == 0) # 31061 show no change across all properties (30.70%)
count(nonsex_specific_snps, sum == 1) # 45645 show one change across all properties (45.12%)
count(nonsex_specific_snps, sum == 2) # 21766 show two changes across all properties (21.9%)
count(nonsex_specific_snps, sum == 3) # 2688 show change across all properties (2.66%)

  # 1. Sex chromosome

count(sexchromo_aa_change, sum == 0) # 1624 show no change across all properties (30.1%)
count(sexchromo_aa_change, sum == 1) # 2436 show one change across all properties (45.2%)
count(sexchromo_aa_change, sum == 2) # 1183 show two changes across all properties (21.52%)
count(sexchromo_aa_change, sum == 3) # 151 show change across all properties (2.8%)

  # 2. Autosomes

count(autosome_aa_change, sum == 0) # 29415 (30.7%)
count(autosome_aa_change, sum == 1) # 43178 (45.1%)
count(autosome_aa_change, sum == 2) # 20569 (21.5%)
count(autosome_aa_change, sum == 3) # 2531 (2.6%)

  # 3. Compared to your missense SNPs that are sex-limited (66):

# 0 changes = 20/66 =
# 1 change = 29/66 =
# 2 changes = 12/66 = 
# 3 changes = 5/66 = 

# Structure = 19/66 = 
# Polarity = 31/66 = 
# Charge = 19/66 =

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

prop.test(snp_rows)


############################################

# 1. No changes

prop.test(c(31061, 19),
          c(101160, 62))
# Output: X-squared = 0.95684, df = 2, p-value = 1

# 2. One change

prop.test(c(45645, 27),
          c(101160, 62))
#Output: X-squared = 0.065367, df = 2, p-value = 0.9035

# 3. Two changes

prop.test(c(21766, 12),
          c(101160, 62))
#Output: X-squared = 0.74913, df = 2, p-value = 0.7953

# 4. All changes

prop.test(c(2688, 4),
          c(101160, 62))
#Output: X-squared = 3.9299, df = 2, p-value = 0.1439



  ### Look at the percent change in the amino acids - individual properties ###

  # 1. Structure

count(autosome_aa_change, structure == 1) # 23196 (24.24%)
count(sexchromo_aa_change, structure == 1) # 1375 (25.5%)
count(nonsex_specific_snps, structure == 1) # 24571 (24.29%)

  # 2. Polarity

count(autosome_aa_change, polarity == 1) # 37423 (39.11%)
count(sexchromo_aa_change, polarity == 1) # 2133 (39.54%)
count(nonsex_specific_snps, polarity == 1) # 39556 (39.10%)

  # 3. Charge

count(autosome_aa_change, charge == 1) # 31290 (32.70%)
count(sexchromo_aa_change, charge == 1) # 1824 (33.82%)
count(nonsex_specific_snps, charge == 1) # 33114 (32.73)


    #### Calculate Stats ###


# 5. Strucutre

prop.test(c(24571, 16),
          c(101160, 62))
# Output: X-squared = 0.95684, df = 2, p-value = 0.8963

# 6. Polarity

prop.test(c(39556, 31),
          c(101160, 62))
#Output: X-squared = 0.065367, df = 2, p-value = 0.1036

# 7. Charge

prop.test(c(33114, 16),
          c(101160, 62))
#Output: X-squared = 0.74913, df = 2, p-value = 0.3045


######################################################################
############# When split up across chromosome types ##################
######################################################################

# Proportion test -  Number of changes:

  # 1. No changes

prop.test(c(29415, 1624, 19),
          c(95693, 5394, 62))
# Output: X-squared = 0.95684, df = 2, p-value = 0.6198

  # 2. One change

prop.test(c(43178, 2436, 27),
          c(95693, 5394, 62))
#Output: X-squared = 0.065367, df = 2, p-value = 0.9678

  # 3. Two changes

prop.test(c(20569, 1183, 12),
          c(95693, 5394, 62))
#Output: X-squared = 0.74913, df = 2, p-value = 0.6876

  # 4. All changes

prop.test(c(2531, 151, 4),
          c(95693, 5394, 62))
#Output: X-squared = 3.9299, df = 2, p-value = 0.1402


# Proportion test -  Type of changes:

  # 1. Structure

prop.test(c(23196, 1375, 16),
          c(95693, 5394, 62))
# Output: X-squared = 4.4209, df = 2, p-value = 0.1097

  # 2. Polarity

prop.test(c(37423, 2133, 31),
          c(95693, 5394, 62))
#Output: X-squared = 3.4818, df = 2, p-value = 0.1754

  # 3. Charge

prop.test(c(31290, 1824, 16),
          c(95693, 5394, 62))
#Output: X-squared = 4.2521, df = 2, p-value = 0.1193




