#SNP density of P. reticulata FINAL

rm(list=ls())
ls() 

# 1. Load libraries & Main File

library(stringr)
library(dplyr)
library(tidyverse)
library(gsubfn)
library(ggplot2)
library(RColorBrewer)

eb <- element_blank()
rstudioapi::writeRStudioPreference("data_viewer_max_columns", 1000L)

setwd("//files.zoology.ubc.ca/ljmfong/flex/new_Ret_SNPs")
getwd()

##########################################################################################
#### Load up the HIGH impact and MODERATE impact files that were sorted with SnpSift ####
##########################################################################################

high_impact_snps <- read.table(file = "filtered_combo_SNPs/formatted_high_impact_SNPs.vcf", sep = "\t", header = TRUE)
moderate_impact_snps <- read.table(file = "filtered_combo_SNPs/formatted_moderate_impact_SNPs.vcf", sep = "\t", header = TRUE)

    # 2. Clean up the descriptions in the file
#First we'll do the high-impact SNPs

unchanged_columns <- data.frame(high_impact_snps[, 1:2])
cleaned_GT <- data.frame(apply(high_impact_snps[, 10:135], 2, function(x) gsub(":.*", "",x)))
final_GT_high_impact <- cbind(unchanged_columns,cleaned_GT)
final_GT_high_impact$IMPACT <- "high"

#Now we'll do the moderate-impact SNPs

unchanged_columns <- data.frame(moderate_impact_snps[, 1:2])
cleaned_GT <- data.frame(apply(moderate_impact_snps[, 10:135], 2, function(x) gsub(":.*", "",x)))
final_GT_moderate_impact <- cbind(unchanged_columns,cleaned_GT)
final_GT_moderate_impact$IMPACT <- "moderate"

#Combine your two dataframes

combo_snps_noMAF <- rbind(final_GT_high_impact, final_GT_moderate_impact)

#Count the number of SNPs and remove min allele freq. < 5% a.k.a 6 individuals

count_SNPs <- data.frame(rowSums(combo_snps_noMAF == "0/1" | combo_snps_noMAF == "0/2" | combo_snps_noMAF == "1/2"))

combo_snps_unfilt <- cbind(combo_snps_noMAF, count_SNPs)
colnames(combo_snps_unfilt)[130] <- "SNP_count"
combo_snps <- combo_snps_unfilt %>%
  filter(SNP_count > 12) #This will filter for SNPs to present in at least 10% of the individuals

    # 3. Read in chromosomes & Organize SNPs by sexes

chr_sizes <- read.table(file = "//files.zoology.ubc.ca/ljmfong/flex/guppy_chrms_len.txt", sep = "\t", header = TRUE)

#In order to plot your SNPs, you'll want to layer your males and females - therefore it's best two make two dataframes for your plot. In this case, you're going to make one for males (high & moderate impact) and one for females.

male_snps <- cbind(combo_snps[, 1:2], combo_snps[, 13:22], combo_snps[, 33:42], combo_snps[, 53:61],
                   combo_snps[, 72:81], combo_snps[, 92:101], combo_snps[, 112:122], combo_snps[, 126:129])

female_snps <- cbind(combo_snps[, 1:12], combo_snps[, 23:32], combo_snps[, 43:52], combo_snps[, 62:71],
                     combo_snps[, 82:91], combo_snps[, 102:111], combo_snps[, 123:125], combo_snps[, 129])
colnames(female_snps)[66] <- "IMPACT"

male_snps$CHROM <- factor(male_snps$CHROM, paste0('LG', 1:23), paste0('LG', 1:23))
female_snps$CHROM <- factor(female_snps$CHROM, paste0('LG', 1:23), paste0('LG', 1:23))
chr_sizes$CHROM <- factor(chr_sizes$CHROM, paste0('LG', 1:23), paste0('LG', 1:23))


    # 4. Re-organize to account for the number of SNPs present
#You will want to count how many individuals that are female vs males) are heterozygous

fem_heter <- data.frame(rowSums(female_snps == "0/1" | female_snps == "0/2" | female_snps == "1/2"))
male_heter <- data.frame(rowSums(male_snps == "0/1" | male_snps == "0/2" | male_snps == "1/2"))

plot_male <- cbind(male_snps[,1:2], male_snps[,66], male_heter)
colnames(plot_male)[3]<- "IMPACT"
colnames(plot_male)[4] <- "hetero_counts"

plot_female <- cbind(female_snps[,1:2], female_snps[,66], fem_heter)
colnames(plot_female)[3]<- "IMPACT"
colnames(plot_female)[4] <- "hetero_counts"

plot_hetero_ratios <- cbind(male_snps[,1:2], male_snps[,66], male_heter, fem_heter)
colnames(plot_hetero_ratios)[3]<- "IMPACT"
colnames(plot_hetero_ratios)[4]<- "Male_Heterozygous_Counts"
colnames(plot_hetero_ratios)[5]<- "Female_Heterozygous_Counts"

MF_ratio <- log2((plot_male$hetero_counts+0.1)/(plot_female$hetero_counts+0.1))
MF_count_ratio <- (plot_male$hetero_counts+0.1)/(plot_female$hetero_counts+0.1)

plot_hetero_ratios$MF_log2 <- MF_ratio
plot_hetero_ratios$MF_count_ratio <- MF_count_ratio
colnames(plot_hetero_ratios)[6] <- "MF_log2"


    # 5. Remove also the errors that are found in the reference genome 
# e.g. all males  and all females are heterozygous at one site

rm_invariants <- plot_hetero_ratios %>%
  filter_all(any_vars(plot_hetero_ratios$MF_log2 != 0.00000000))

rm_invariants_uniq <- rm_invariants %>%
  distinct(POS, .keep_all = TRUE)

autosomes <- subset(rm_invariants_uniq, CHROM != "LG12")
only_LG12 <- subset(rm_invariants_uniq, CHROM == "LG12")

range(rm_invariants$MF_log2)
#[1] -7.238405  8.535275

#Clean up a few files:
rm(unchanged_columns)
rm(cleaned_GT)
rm(combo_snps_noMAF)
rm(count_SNPs)
rm(combo_snps_unfilt)
rm(male_snps)
rm(female_snps)
rm(plot_female)
rm(plot_male)
rm(fem_heter)
rm(male_heter)

    # 6. Base Plots:

chr_sizes$start <- 0

base_wouter <- ggplot() +
  geom_rect(data = chr_sizes, aes(xmin = 0, xmax = LEN, ymin = 0, ymax = 25), col = NA, fill = 'grey95', orientation = 'y') +
  geom_rect(data = subset(chr_sizes, CHROM == "LG12"), aes(xmin = 20e6, xmax = 26e6, ymin = 0, ymax = 25), col = NA, fill = 'mediumpurple3', orientation = 'y', alpha = 0.3) +
  geom_segment(aes(x = 0, xend = LEN, y = 0, yend = 0), chr_sizes, col = "black") +
  scale_x_continuous(expand = c(0,0), labels = ~.x / 1e6, name = 'Chromosomal position (Mb)') +
  scale_y_continuous(lim = c(0, 26), breaks = c()) +
  facet_grid(rows = vars(CHROM), switch = 'y') +
  theme_classic() +
  theme(
    legend.position = 'top', strip.placement = 'outside',
    strip.background = eb, strip.text.y.left = element_text(angle = 0,size = 18),
    axis.text.x = element_text(size = 22), axis.text.y = element_text(size = 18), axis.title = element_text(size = 27)
  ) + ylab("Chromosomes") +  xlab("Chromosome (Mb)") +
  labs(alpha = "No. of Individuals with SNP")
                               
chrm_12 <- subset(chr_sizes, CHROM == "LG12")

    # 7. Organize sex-limited SNPs:

male_limited_snps <- rm_invariants %>%
  filter_all(any_vars(rm_invariants$Female_Heterozygous_Counts == 0)) 

female_limited_snps <- rm_invariants %>%
  filter_all(any_vars(rm_invariants$Male_Heterozygous_Counts == 0)) 

#You may get repeats due to how SNPeff categorizes impacts, so ensure that your get only unique values (i.e. unique positions)

male_limited_snps_uniq <- male_limited_snps %>%
  distinct(POS, .keep_all = TRUE) #I have 90 male-limited SNPs
write.table(male_limited_snps_uniq, file = "all_male_limited_SNPs.txt", sep = "\t", row.names = TRUE, quote = F)

female_limited_snps_uniq <- female_limited_snps %>%
  distinct(POS, .keep_all = TRUE) # I have 1 female-limited genes

sex_limited_all_snps_uniq <- rbind(male_limited_snps_uniq, female_limited_snps_uniq)
write.table(sex_limited_all_snps_uniq, file = "All_sex_limited_SNPs.txt", sep = "\t", row.names = TRUE, quote = F)

male_limited_LG12 <- subset(male_limited_snps_uniq, CHROM == "LG12")
female_limited_LG12 <- subset(female_limited_snps_uniq, CHROM == "LG12")


    # 8. Plotting

base_wouter +
  geom_point(
    aes(x = POS, y = MF_log2, alpha = Female_Heterozygous_Counts), col = 'red',
    data = female_limited_snps_uniq, size = 3
  ) +
  geom_point(
    aes(x = POS, y = MF_log2, alpha = Male_Heterozygous_Counts), col = 'blue',
    data = male_limited_snps_uniq, size = 3
  )

#### Filter out the impacted genes to check the type of warnings ####

to_filter_high_impact <- high_impact_snps %>%
  semi_join(sex_limited_all_snps_uniq, by = c("CHROM", "POS"))
to_filter_mod_impact <- moderate_impact_snps %>%
  semi_join(sex_limited_all_snps_uniq, by = c("CHROM", "POS"))

to_filter_SNPs <- rbind(to_filter_high_impact, to_filter_mod_impact)
write.table(to_filter_SNPs, "unfiltered_for_WARNINGS_snps.txt", sep = "\t", row.names = T, quote = F)
# If there are any WARNINGs that suggest there is an issue with the SNP prediction, remove it from downstream analyses

########### Lollipop Plots: #############

#Bar Plot + all chrms plot:
plot_LG12_lollipop <- ggplot(only_LG12, aes(x= POS/1000000)) +
  scale_x_continuous(lim = c(0, 27000000), breaks = c(0, 5e6, 10e6, 15e6, 20e6, 25e6), 
                     labels = ~.x / 1e6, name = 'Chromosome 12 (Mb)') +
  scale_y_continuous(expand = c(0,0), lim = c(0, 26)) + 
  theme_classic() +
  theme(
    legend.position = 'top', strip.placement = 'outside',
    strip.background = eb, strip.text.y.left = element_text(angle = 0, size = 16), 
    axis.line = element_line(colour = "black", size = 1.2), axis.ticks = element_line(size = 1.5),
    axis.text = element_text(size = 18), axis.title = element_text(size = 20)
  ) + ylab("No. of Individuals with SNP") + xlab("Chromosome 12 (Mb)") + 
  labs(values = c("Male", "Female"))

plot_LG12_lollipop +
  annotate("rect", xmin = 20e6, xmax = 26e6, ymin = 0, ymax = 25, fill = 'mediumpurple3', alpha = 0.3) +
  geom_segment(data = male_limited_LG12,
               aes(x = POS, y = 0, yend = (Male_Heterozygous_Counts-0.5)), size = 1, col = "blue4") +
  geom_point(data = male_limited_LG12, aes(x = POS, y = Male_Heterozygous_Counts),
             size = 4, col = "navy", fill = alpha("cornflowerblue", 0.8), shape = 23, stroke = 1.5) +
  geom_segment(data = female_limited_LG12, 
               aes(x = POS, y = 0, yend = Female_Heterozygous_Counts-0.5), size = 1, col = "firebrick4") +
  geom_point(data = female_limited_LG12, aes(x = POS, y = Female_Heterozygous_Counts),
             size = 4, col = "firebrick4", fill = alpha("firebrick3", 0.8), shape = 21, stroke = 1.5)

#### For all chromosomes ####

base_wouter +
  geom_linerange(data = male_limited_snps_uniq, 
                 aes(x = POS, ymin = 0, ymax = (Male_Heterozygous_Counts-0.5)), size = 1, col = "blue4") +
  geom_point(data = male_limited_snps_uniq, aes(x = POS, y = Male_Heterozygous_Counts),
             size = 3, col = "navy", fill = alpha("cornflowerblue", 0.8), shape = 23, stroke = 1.2) +
  geom_linerange(data = female_limited_snps_uniq, 
                 aes(x = POS, ymin = 0, ymax = (Female_Heterozygous_Counts-0.5)), size = 1, col = "red4") +
  geom_point(data = female_limited_snps_uniq, aes(x = POS, y = Female_Heterozygous_Counts),
             size = 3, col = "firebrick4", fill = alpha("firebrick3", 0.8), shape = 23, stroke = 1.2)

##########################################################################################
################ Check for dupes and remove those genes and their SNPs! ##################
##########################################################################################

#See threshold_for_MF_depth_final.R to identify genes that are likely duplicated genes

# From the autosomes:
# LG6: ENSPREG00000007572; LG7: ENSPREG00000003921; LG7: ENSPREG00000012441; LG15: ENSPREG00000005754;
# From the sex chromosomes:
# ENSPREG00000019378; ENSPREG00000019774; ENSPREG00000006331; ENSPREG00000006362; ST8SIA3; CCNG2

male_limited_snps_no_dupe <- male_limited_snps_uniq %>%
  filter(CHROM != "LG6" & CHROM != "LG7" & CHROM != "LG15" & POS != 20712005 & POS != 21196665 & POS != 25783963 & POS != 25833552 & POS != 25836614 & POS != 25840579 & POS != 21214735 & POS != 21348667)

male_limited_snps_no_dupe_LG12 <- subset(male_limited_snps_no_dupe, CHROM == "LG12")

plot_LG12_lollipop +
  annotate("rect", xmin = 20e6, xmax = 26e6, ymin = 0, ymax = 26, fill = 'mediumpurple3', alpha = 0.3) +
  geom_segment(data = male_limited_snps_no_dupe_LG12,
               aes(x = POS, y = 0, yend = (Male_Heterozygous_Counts-0.5)), size = 1, col = "blue4") +
  geom_point(data = male_limited_snps_no_dupe_LG12, aes(x = POS, y = Male_Heterozygous_Counts),
             size = 5, col = "navy", fill = alpha("cornflowerblue", 0.8), shape = 21, stroke = 1.5) +
  geom_segment(data = female_limited_LG12, 
               aes(x = POS, y = 0, yend = Female_Heterozygous_Counts), size = 1) +
  geom_point(data = female_limited_LG12, aes(x = POS, y = Female_Heterozygous_Counts),
             size = 5, col = "firebrick4", fill = alpha("firebrick3", 0.8), shape = 21, stroke = 1.5)

base_wouter +
  geom_linerange(data = male_limited_snps_no_dupe, 
                 aes(x = POS, ymin = 0, ymax = (Male_Heterozygous_Counts-3)), size = 1, col = "blue4") +
  geom_point(data = male_limited_snps_no_dupe, aes(x = POS, y = Male_Heterozygous_Counts),
             size = 4, col = "navy", fill = alpha("cornflowerblue", 0.8), shape = 21, stroke = 1.2) +
  geom_linerange(data = female_limited_snps_uniq, 
                 aes(x = POS, ymin = 0, ymax = (Female_Heterozygous_Counts-3)), size = 1, col = "red4") +
  geom_point(data = female_limited_snps_uniq, aes(x = POS, y = Female_Heterozygous_Counts),
             size = 4, col = "firebrick4", fill = alpha("firebrick3", 0.8), shape = 21, stroke = 1.2)


### When you know which positions have your high impact SNPs after filtering - pull out the sequence and position:

# On LG12 (moderate): 5652422, 5705274, 5705331,
# 20712005, 20931822, 20936465, 20936621, 20938276, 20942935, 21072080, 21214735, 21285375,
# 21285399, 21328421, 21340862, 21344523, 21344563, 21345060, 21348667, 21359254, 21364902, 21365941, 21372888,
# 21385522, 21385618, 21385891, 21387147, 21387430, 21387661, 21391427, 21391639, 21438860, 21443385, 21466508,
# 21466550, 21468271, 21468532, 21486526, 21513934, 21514479, 21514515, 21515845, 21888711, 24831338, 24923220,
# 24964412, 24964761, 24996712, 24998564, 25141492, 25141630, 25141935, 25141951, 25261456, 25261461, 25267942,
# 25275647, 25278301, 25280384, 25280689, 25281082, 25281125, 25281231, 25281902, 25281946, 25319803, 25319966,
# 25321077, 25833552, 25836614, 25840579, 25861233, 25863124

# On LG12 (high): 21050147, 21343808, 24930241

mod_list <- c(5652422, 5705274, 5705331, 20712005, 20931822, 20936465, 20936621, 20938276, 20942935, 21072080, 21214735, 21285375,
              21285399, 21328421, 21340862, 21344523, 21344563, 21345060, 21348667, 21359254, 21364902, 21365941, 21372888,
              21385522, 21385618, 21385891, 21387147, 21387430, 21387661, 21391427, 21391639, 21438860, 21443385, 21466508,
              21466550, 21468271, 21468532, 21486526, 21513934, 21514479, 21514515, 21515845, 21888711, 24831338, 24923220,
              24964412, 24964761, 24996712, 24998564, 25141492, 25141630, 25141935, 25141951, 25261456, 25261461, 25267942,
              25275647, 25278301, 25280384, 25280689, 25281082, 25281125, 25281231, 25281902, 25281946, 25319803, 25319966,
              25321077, 25833552, 25836614, 25840579, 25861233, 25863124)
high_list <- c(21050147, 21343808, 24930241)

snp_change_high <- filter(high_impact_snps, grepl(paste(high_list, collapse = "|"), POS))
snp_change_mod <- filter(moderate_impact_snps, grepl(paste(mod_list, collapse = "|"), POS))
combo_snp_change <- rbind(snp_change_high, snp_change_mod)

write.table(combo_snp_change, "snp_change_LG12.txt", sep = "\t", quote = F)



