#### Plotting the ase between males and females from scRNA ####

rm(list=ls())
ls() 

setwd("//files.zoology.ubc.ca/ljmfong/flex/cellSNP_scAlleleCount/R_results")
getwd()

#Ensure you're on the latest R version

#Load in libraries

library(ggplot2)
library(dplyr)
library(patchwork)
library(rstatix)
library(purrr)
library(broom)
library(tidyverse)
library(patchwork)
library(ggpubr)


##############################################################################
########### Pseudobulk and look at count distribution of tissues #############
##############################################################################

gonad_minfilt <- read.table(file = "snpfilt_gonad_filtered10.txt", header = T)
skin_minfilt <- read.table(file = "snpfilt_skin_filtered10.txt", header = T)
liver_minfilt <- read.table(file = "snpfilt_liver_filtered10.txt", header = T)
heart_minfilt <- read.table(file = "snpfilt_heart_filtered10.txt", header = T)

gonad_minfilt$tissue <- "gonad"
skin_minfilt$tissue <- "skin"
liver_minfilt$tissue <- "liver"
heart_minfilt$tissue <- "heart"

gonad_minfilt$sex <- ifelse(grepl("^F", gonad_minfilt$sample_id), "F", "M")
skin_minfilt$sex <- ifelse(grepl("^F", skin_minfilt$sample_id), "F", "M")
liver_minfilt$sex <- ifelse(grepl("^F", liver_minfilt$sample_id), "F", "M")
heart_minfilt$sex <- ifelse(grepl("^F", heart_minfilt$sample_id), "F", "M")

somatic_all <- rbind(skin_minfilt, liver_minfilt, heart_minfilt)
all_tissue <- rbind(somatic_all, gonad_minfilt)


#################################################################
######### Look at a cell-by-cell rather than average ############
#################################################################

# Attach gene information and chromosome type 

gene_info <- read.table("\\\\files.zoology.ubc.ca/ljmfong/flex/new_Ret_SNPs/gene_boundary_w_gene_names.bed", header = T, sep = "\t")
gene_list <- c("unc5da", "ENSPREG00000019478", "ENSPREG00000019540", "ENSPREG00000019560", "ENSPREG00000019673",
               "ENSPREG00000019838", "WDR7", "ENSPREG00000019881", "ccni", "ENSPREG00000019970", "ENSPREG00000019988", 
               "ENSPREG00000020020", "ENSPREG00000020032", "ENSPREG00000020049", "ENSPREG00000020126", "ENSPREG00000020131", 
               "zfand5b", "ccnb1", "CACNA1B", "ENSPREG00000001200", "c7b", "ENSPREG00000001446", "ENSPREG00000001498",
               "ENSPREG00000005507", "ENSPREG00000005501", "ENSPREG00000005616", "PIGO", "ENSPREG00000006387") 

non_recombo_gene_list <- gene_info %>%
  filter(CHROM == "LG12", START >= 19992643, START <= 25999361) %>% pull(GENE_NAME)
#Only about 174 genes are in the non-recombining region

add_chrom_column <- function(dfA, dfB) {
  matched_idx <- match(dfA$geneid, dfB$GENE_NAME)
  dfA$chrom <- dfB$CHROM[matched_idx]
  return(dfA)
}

all_tissue <- add_chrom_column(all_tissue, gene_info)
all_tissue$chrom_type <- ifelse(all_tissue$chrom == "LG12", "Sex Chromosome", "Autosome")
all_tissue$chrom_type[!all_tissue$geneid %in% non_recombo_gene_list] <- "Autosome"
all_tissue$chrom_type[all_tissue$geneid %in% gene_list] <- "Impacted Gene"

#################################################################################
################ Filter for ref strand bias in different samples ################
#################################################################################

somatic_all <- all_tissue %>% filter(tissue != "gonad")
gonad_all <- all_tissue %>% filter(tissue == "gonad")

somatic_allele_counts <- somatic_all
gonad_allele_counts <- gonad_all


#################################################################
########### Organize tissues, then plot and stats ###############
#################################################################

# GONAD #

gonad_ase <- gonad_allele_counts
gonad_ase$total_counts <- gonad_ase$ref_count + gonad_ase$alt_count
gonad_ase$mar <- pmax(gonad_ase$ref_count, gonad_ase$alt_count)/(gonad_ase$total_counts)
gonad_ase <- na.omit(gonad_ase)

gonad_samp_filt <- gonad_ase %>% group_by(geneid) %>%
  filter(sum(sex == "F") >= 3, sum(sex == "M") >= 3) %>% ungroup()

gonad_autosome <- subset(gonad_samp_filt, chrom_type == "Autosome")
gonad_sexchromo <- subset(gonad_samp_filt, chrom_type == "Sex Chromosome")

gonad_autosome %>% filter(sex == "F") %>% pull(mar) %>% median() # 0.5971564
gonad_autosome %>% filter(sex == "M") %>% pull(mar) %>% median() # 0.5862069
gonad_sexchromo %>% filter(sex == "F") %>% pull(mar) %>% median() # 0.5893997
gonad_sexchromo %>% filter(sex == "M") %>% pull(mar) %>% median() # 0.5833333

plot_gonad_A <- ggplot(gonad_autosome, aes(x = mar)) + geom_density(aes(fill = sex, col = sex), alpha = 0.5, adjust = 1.5, linewidth = 1.2) + 
  scale_fill_manual(values = c("firebrick3", "cornflowerblue")) +  labs(fill = "Sex", colour = "Sex") +
  scale_colour_manual(values = c("firebrick3", "cornflowerblue")) +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) + scale_y_continuous(expand = c(0,0), lim = c(0, 7)) +
  ylab("Density") + xlab("") + ggtitle("") +
  geom_vline(aes(xintercept=0.5971564),color="firebrick4",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5862069),color="navy",linetype="dashed",size=0.8) +
  theme_classic() + theme( legend.position = "top",
    axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18),
    plot.title = element_text(size = 20, hjust = 0.5)
  )
plot_gonad_A

plot_gonad_SC <- ggplot(gonad_sexchromo, aes(x = mar)) + geom_density(aes(fill = sex, col = sex), alpha = 0.5, adjust = 1.5, linewidth = 1.2) + 
  scale_fill_manual(values = c("firebrick3", "cornflowerblue")) +  labs(fill = "Sex", colour = "Sex") +
  scale_colour_manual(values = c("firebrick3", "cornflowerblue")) +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) + scale_y_continuous(expand = c(0,0), lim = c(0, 7)) +
  ylab("Density") + xlab("") + ggtitle("") +
  geom_vline(aes(xintercept=0.5893997),color="firebrick4",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5833333),color="navy",linetype="dashed",size=0.8) +
  theme_classic() + theme( legend.position = "top",
    axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),  axis.text.y = element_text(size = 18),
    plot.title = element_text(size = 20, hjust = 0.5)
  )
plot_gonad_SC

gonad_two_plots <- (plot_gonad_A | plot_gonad_SC) + plot_layout(guides = "collect")
gonad_two_plots

chisq_test(table(gonad_autosome$sex, gonad_autosome$mar))
# 0.0138
chisq_test(table(gonad_sexchromo$sex, gonad_sexchromo$mar))
#  0.414


# SOMATIC #

somatic_ase <- somatic_allele_counts
somatic_ase$total_counts <- somatic_ase$ref_count + somatic_ase$alt_count
somatic_ase$mar <- pmax(somatic_ase$ref_count, somatic_ase$alt_count)/(somatic_ase$total_counts)
somatic_ase <- na.omit(somatic_ase)

somatic_samp_filt <- somatic_ase %>% group_by(geneid) %>%
  filter(sum(sex == "F") >= 3, sum(sex == "M") >= 3) %>% ungroup()

somatic_autosome <- subset(somatic_samp_filt, chrom_type == "Autosome")
somatic_sexchromo <- subset(somatic_samp_filt, chrom_type == "Sex Chromosome")

somatic_autosome %>% filter(sex == "F") %>% pull(mar) %>% median() # 0.5862069
somatic_autosome %>% filter(sex == "M") %>% pull(mar) %>% median() # 0.6078431
somatic_sexchromo %>% filter(sex == "F") %>% pull(mar) %>% median() # 0.575
somatic_sexchromo %>% filter(sex == "M") %>% pull(mar) %>% median() # 0.617

plot_somatic_A <- ggplot(somatic_autosome, aes(x = mar)) + geom_density(aes(fill = sex, col = sex), alpha = 0.5, adjust = 1.5, linewidth = 1.2) + 
  scale_fill_manual(values = c("firebrick3", "cornflowerblue")) + labs(fill = "Sex", colour = "Sex") +
  scale_colour_manual(values = c("firebrick3", "cornflowerblue")) +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +  scale_y_continuous(expand = c(0,0), lim = c(0, 7)) +
  ylab("Density") + xlab("Major Allele Ratio") + ggtitle("") +
  geom_vline(aes(xintercept=0.5862069),color="firebrick4",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6078431),color="navy",linetype="dashed",size=0.8) +
  theme_classic() + theme( legend.position = "top",
    axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 18),  axis.text.y = element_text(size = 18),
    plot.title = element_text(size = 20, hjust = 0.5)
  )
plot_somatic_A

plot_somatic_SC <- ggplot(somatic_sexchromo, aes(x = mar)) + geom_density(aes(fill = sex, col = sex), alpha = 0.5, adjust = 1.5, linewidth = 1.2) + 
  scale_fill_manual(values = c("firebrick3", "cornflowerblue")) + labs(fill = "Sex", colour = "Sex") +
  scale_colour_manual(values = c("firebrick3", "cornflowerblue")) +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +  scale_y_continuous(expand = c(0,0), lim = c(0, 7)) +
  ylab("Density") + xlab("Major Allele Ratio") + ggtitle("") +
  geom_vline(aes(xintercept=0.575),color="firebrick4",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.617),color="navy",linetype="dashed",size=0.8) +
  theme_classic() + theme( legend.position = "top",
    axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 18),  axis.text.y = element_text(size = 18),
    plot.title = element_text(size = 20, hjust = 0.5)
  )
plot_somatic_SC

somatic_two_plots <- (plot_somatic_A | plot_somatic_SC) + plot_layout(guides = "collect")
somatic_two_plots

chisq_test(table(somatic_autosome$sex, somatic_autosome$mar))
# 0.253
chisq_test(table(somatic_sexchromo$sex, somatic_sexchromo$mar))
# 0.868 

combo_tissue_plots <- (plot_gonad_A | plot_gonad_SC) / (plot_somatic_A | plot_somatic_SC) + plot_layout(guides = "collect")
combo_tissue_plots


##################################################################
########## Plot each tissue type by sex and chromo type ##########
##################################################################

somatic_ase <- somatic_allele_counts
somatic_ase$chrom_type <- factor(somatic_ase$chrom_type, levels = c("Autosome", "Sex Chromosome", "Impacted Gene"))
somatic_ase$total_counts <- somatic_ase$ref_count + somatic_ase$alt_count
somatic_ase$mar <- pmax(somatic_ase$ref_count, somatic_ase$alt_count)/(somatic_ase$total_counts)
somatic_ase <- na.omit(somatic_ase)

gonad_ase <- gonad_allele_counts
gonad_ase$chrom_type <- factor(gonad_ase$chrom_type, levels = c("Autosome", "Sex Chromosome", "Impacted Gene"))
gonad_ase$total_counts <- gonad_ase$ref_count + gonad_ase$alt_count
gonad_ase$mar <- pmax(gonad_ase$ref_count, gonad_ase$alt_count)/(gonad_ase$total_counts)
gonad_ase <- na.omit(gonad_ase)

male_somatic_allele_counts <- subset(somatic_ase, sex == "M")
male_somatic_allele_counts <- na.omit(male_somatic_allele_counts)
male_gonad_allele_counts <- subset(gonad_ase, sex == "M")
male_gonad_allele_counts <- na.omit(male_gonad_allele_counts)
female_somatic_allele_counts <- subset(somatic_ase, sex == "F")
female_somatic_allele_counts <- na.omit(female_somatic_allele_counts)
female_gonad_allele_counts <- subset(gonad_ase, sex == "F")
female_gonad_allele_counts <- na.omit(female_gonad_allele_counts)

# MALE #

clean_male_somatic_allele_counts <- subset(male_somatic_allele_counts, chrom_type != "Impacted Gene")
clean_male_somatic_allele_counts <- droplevels(clean_male_somatic_allele_counts, "Impacted Gene")
clean_male_somatic_allele_counts <- clean_male_somatic_allele_counts %>% group_by(geneid) %>%
  filter(sum(sex == "M") >= 3) %>% ungroup()

clean_male_somatic_allele_counts %>% filter(chrom_type == "Autosome") %>% pull(mar) %>% median() # 0.6071429
clean_male_somatic_allele_counts %>% filter(chrom_type == "Sex Chromosome") %>% pull(mar) %>% median() # 0.617

clean_plot_somatic_male <- ggplot(clean_male_somatic_allele_counts, aes(x = mar)) + geom_density(aes(fill = chrom_type, col = chrom_type), alpha = 0.6, adjust = 1.5, linewidth = 1.2) + 
  scale_fill_manual(values = c("#857C7B", "mediumpurple3")) + labs(fill = "SNP Type", col = "SNP Type") +
  scale_colour_manual(values = c("#857C7B", "mediumpurple3")) + 
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +  scale_y_continuous(expand = c(0,0), lim = c(0, 7)) +
  ylab("") + xlab("") + ggtitle("") +
  geom_vline(aes(xintercept=0.6071429),color="grey40",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.617),color="darkorchid3",linetype="dashed",size=0.8) +
  theme_classic() + theme(
    axis.title.x = element_text(size = 25),  axis.title.y = element_text(size = 25), legend.position = "top",
    axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
    plot.title = element_text(size = 20, hjust = 0.5),
    plot.margin = margin(t = 20, r = 50, b = 50, l = 20)
  )
clean_plot_somatic_male

clean_male_gonad_allele_counts <- subset(male_gonad_allele_counts, chrom_type != "Impacted Gene")
clean_male_gonad_allele_counts<- droplevels(clean_male_gonad_allele_counts, "Impacted Gene")
clean_male_gonad_allele_counts <- clean_male_gonad_allele_counts %>% group_by(geneid) %>%
  filter(sum(sex == "M") >= 3) %>% ungroup()

clean_male_gonad_allele_counts %>% filter(chrom_type == "Autosome") %>% pull(mar) %>% median() # 0.5862069
clean_male_gonad_allele_counts %>% filter(chrom_type == "Sex Chromosome") %>% pull(mar) %>% median() # 0.5777778

clean_plot_gonad_male <- ggplot(clean_male_gonad_allele_counts, aes(x = mar)) + geom_density(aes(fill = chrom_type, col = chrom_type), alpha = 0.5, adjust = 1.5, linewidth = 1.2) + 
  scale_fill_manual(values = c("#857C7B", "mediumpurple3")) + labs(fill = "SNP Type", col = "SNP Type") +
  scale_colour_manual(values = c("#857C7B", "mediumpurple3")) + 
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +  scale_y_continuous(expand = c(0,0), lim = c(0, 7)) +
  ylab("") + xlab("") + ggtitle("") +
  geom_vline(aes(xintercept=0.5862069),color="grey40",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5777778),color="darkorchid3",linetype="dashed",size=0.8) +
  theme_classic() + theme(
    axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25), legend.position = "top",
    axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
    plot.title = element_text(size = 20, hjust = 0.5),
    plot.margin = margin(t = 20, r = 50, b = 50, l = 20)
  )
clean_plot_gonad_male

chisq_test(table(clean_male_somatic_allele_counts$chrom_type, clean_male_somatic_allele_counts$mar))
# 1.00
chisq_test(table(clean_male_gonad_allele_counts$chrom_type, clean_male_gonad_allele_counts$mar))
# 1.00

# FEMALE #

clean_female_somatic_allele_counts <- subset(female_somatic_allele_counts, chrom_type != "Impacted Gene")
clean_female_somatic_allele_counts<- droplevels(clean_female_somatic_allele_counts, "Impacted Gene")
clean_female_somatic_allele_counts <- clean_female_somatic_allele_counts %>% group_by(geneid) %>%
  filter(sum(sex == "F") >= 3) %>% ungroup()

clean_female_somatic_allele_counts %>% filter(chrom_type == "Autosome") %>% pull(mar) %>% median() # 0.5862069
clean_female_somatic_allele_counts %>% filter(chrom_type == "Sex Chromosome") %>% pull(mar) %>% median() # 0.5714286

clean_plot_somatic_female <- ggplot(clean_female_somatic_allele_counts, aes(x = mar)) + geom_density(aes(fill = chrom_type, col = chrom_type), alpha = 0.5, adjust = 1.5, linewidth = 1.2) + 
  scale_fill_manual(values = c("#857C7B", "mediumpurple3")) + labs(fill = "SNP Type", col = "SNP Type") +
  scale_colour_manual(values = c("#857C7B", "mediumpurple3")) + 
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +  scale_y_continuous(expand = c(0,0), lim = c(0, 7)) +
  ylab("") + xlab("") + ggtitle("") +
  geom_vline(aes(xintercept=0.5862069),color="grey40",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5714286),color="darkorchid3",linetype="dashed",size=0.8) +
  theme_classic() + theme(
    axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25), legend.position = "none",
    axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
    plot.title = element_text(size = 20, hjust = 0.5),
    plot.margin = margin(t = 20, r = 50, b = 50, l = 20)
  )
clean_plot_somatic_female

clean_female_gonad_allele_counts <- subset(female_gonad_allele_counts, chrom_type != "Impacted Gene")
clean_female_gonad_allele_counts<- droplevels(clean_female_gonad_allele_counts, "Impacted Gene")
clean_female_gonad_allele_counts <- clean_female_gonad_allele_counts %>% group_by(geneid) %>%
  filter(sum(sex == "F") >= 3) %>% ungroup()

clean_female_gonad_allele_counts %>% filter(chrom_type == "Autosome") %>% pull(mar) %>% median() # 0.5957447
clean_female_gonad_allele_counts %>% filter(chrom_type == "Sex Chromosome") %>% pull(mar) %>% median() # 0.6060606

clean_plot_gonad_female <- ggplot(clean_female_gonad_allele_counts, aes(x = mar)) + geom_density(aes(fill = chrom_type, col = chrom_type), alpha = 0.5, adjust = 1.5, linewidth = 1.2) + 
  scale_fill_manual(values = c("#857C7B", "mediumpurple3")) + labs(fill = "SNP Type", col = "SNP Type") +
  scale_colour_manual(values = c("#857C7B", "mediumpurple3")) + 
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +  scale_y_continuous(expand = c(0,0), lim = c(0, 7)) +
  ylab("") + xlab("") + ggtitle("") +
  geom_vline(aes(xintercept=0.5957447),color="grey40",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6060606),color="darkorchid3",linetype="dashed",size=0.8) +
  theme_classic() + theme(
    axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25), legend.position = "none",
    axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
    plot.title = element_text(size = 20, hjust = 0.5),
    plot.margin = margin(t = 20, r = 50, b = 50, l = 20)
  )
clean_plot_gonad_female

chisq_test(table(clean_female_somatic_allele_counts$chrom_type, clean_female_somatic_allele_counts$mar), correct = T)
# 1 
chisq_test(table(clean_female_gonad_allele_counts$chrom_type, clean_female_gonad_allele_counts$mar), correct = T)
# 0.108 

snp_type_plots <- (clean_plot_somatic_female | clean_plot_gonad_female) / (clean_plot_somatic_male | clean_plot_gonad_male) +
  plot_layout(guides = "collect")
snp_type_plots


##############################################################################
################### Plot all impacted genes male alleles: ####################
##############################################################################

#### Do this first for Somatic Tissues ####

M_somatic_CACNA1B_df <- somatic_allele_counts %>% filter(geneid == "CACNA1B" & sex == "M")
M_somatic_CACNA1B_df$total_count <- M_somatic_CACNA1B_df$ref_count + M_somatic_CACNA1B_df$alt_count
M_somatic_CACNA1B_df_long <- M_somatic_CACNA1B_df %>%
  pivot_longer(cols = c(ref_count, alt_count), names_to = "count_type", values_to = "count") %>%
  mutate(count_type = case_when(count_type == "ref_count" ~ "Reference", count_type == "alt_count" ~ "Alternate"))
M_somatic_CACNA1B_df_long <- M_somatic_CACNA1B_df_long [order(M_somatic_CACNA1B_df_long$count, decreasing = F), ]

M_somatic_ENSPREG00000001200_df <- somatic_allele_counts %>% filter(geneid == "ENSPREG00000001200" & sex == "M")
M_somatic_ENSPREG00000001200_df$total_count <- M_somatic_ENSPREG00000001200_df$ref_count + M_somatic_ENSPREG00000001200_df$alt_count
M_somatic_ENSPREG00000001200_df_long <- M_somatic_ENSPREG00000001200_df %>%
  pivot_longer(cols = c(ref_count, alt_count), names_to = "count_type", values_to = "count") %>%
  mutate(count_type = case_when(count_type == "ref_count" ~ "Reference", count_type == "alt_count" ~ "Alternate"))
M_somatic_ENSPREG00000001200_df_long <- M_somatic_ENSPREG00000001200_df_long [order(M_somatic_ENSPREG00000001200_df_long$count, decreasing = F), ]

M_somatic_ENSPREG00000005616_df <- somatic_allele_counts %>% filter(geneid == "ENSPREG00000005616" & sex == "M")
M_somatic_ENSPREG00000005616_df$total_count <- M_somatic_ENSPREG00000005616_df$ref_count + M_somatic_ENSPREG00000005616_df$alt_count
M_somatic_ENSPREG00000005616_df_long <- M_somatic_ENSPREG00000005616_df %>%
  pivot_longer(cols = c(ref_count, alt_count), names_to = "count_type", values_to = "count") %>%
  mutate(count_type = case_when(count_type == "ref_count" ~ "Reference", count_type == "alt_count" ~ "Alternate"))
M_somatic_ENSPREG00000005616_df_long <- M_somatic_ENSPREG00000005616_df_long [order(M_somatic_ENSPREG00000005616_df_long$count, decreasing = F), ]

M_somatic_ENSPREG00000019478_df <- somatic_allele_counts %>% filter(geneid == "ENSPREG00000019478" & sex == "M")
M_somatic_ENSPREG00000019478_df$total_count <- M_somatic_ENSPREG00000019478_df$ref_count + M_somatic_ENSPREG00000019478_df$alt_count
M_somatic_ENSPREG00000019478_df_long <- M_somatic_ENSPREG00000019478_df %>%
  pivot_longer(cols = c(ref_count, alt_count), names_to = "count_type", values_to = "count") %>%
  mutate(count_type = case_when(count_type == "ref_count" ~ "Reference", count_type == "alt_count" ~ "Alternate"))
M_somatic_ENSPREG00000019478_df_long <- M_somatic_ENSPREG00000019478_df_long [order(M_somatic_ENSPREG00000019478_df_long$count, decreasing = F), ]

M_somatic_unc5da_df <- somatic_allele_counts %>% filter(geneid == "unc5da" & sex == "M")
M_somatic_unc5da_df$total_count <- M_somatic_unc5da_df$ref_count + M_somatic_unc5da_df$alt_count
M_somatic_unc5da_df_long <- M_somatic_unc5da_df %>%
  pivot_longer(cols = c(ref_count, alt_count), names_to = "count_type", values_to = "count") %>%
  mutate(count_type = case_when(count_type == "ref_count" ~ "Reference", count_type == "alt_count" ~ "Alternate"))
M_somatic_unc5da_df_long <- M_somatic_unc5da_df_long [order(M_somatic_unc5da_df_long$count, decreasing = F), ]

wilcox.test(M_somatic_CACNA1B_df_long$count ~ M_somatic_CACNA1B_df_long$count_type, p.adjust.methods = "bonferroni") # 0.07083
wilcox.test(M_somatic_ENSPREG00000001200_df_long$count ~ M_somatic_ENSPREG00000001200_df_long$count_type, p.adjust.methods = "bonferroni") # 0.1266
wilcox.test(M_somatic_ENSPREG00000005616_df_long$count ~ M_somatic_ENSPREG00000005616_df_long$count_type, p.adjust.methods = "bonferroni") # 0.5597
wilcox.test(M_somatic_ENSPREG00000019478_df_long$count ~ M_somatic_ENSPREG00000019478_df_long$count_type, p.adjust.methods = "bonferroni") # 0.2034
wilcox.test(M_somatic_unc5da_df_long$count ~ M_somatic_unc5da_df_long$count_type, p.adjust.methods = "bonferroni") # 0.002509

#### For Female ####

F_somatic_CACNA1B_df <- somatic_allele_counts %>% filter(geneid == "CACNA1B" & sex == "F")
F_somatic_CACNA1B_df_long <- F_somatic_CACNA1B_df %>%
  pivot_longer(cols = c(ref_count, alt_count), names_to = "count_type", values_to = "count") %>%
  mutate(count_type = case_when(count_type == "ref_count" ~ "Reference", count_type == "alt_count" ~ "Alternate"))
F_somatic_CACNA1B_df_long <- F_somatic_CACNA1B_df_long [order(F_somatic_CACNA1B_df_long$count, decreasing = F), ]

F_somatic_ENSPREG00000001200_df <- somatic_allele_counts %>% filter(geneid == "ENSPREG00000001200" & sex == "F")
F_somatic_ENSPREG00000001200_df_long <- F_somatic_ENSPREG00000001200_df %>%
  pivot_longer(cols = c(ref_count, alt_count), names_to = "count_type", values_to = "count") %>%
  mutate(count_type = case_when(count_type == "ref_count" ~ "Reference", count_type == "alt_count" ~ "Alternate"))
F_somatic_ENSPREG00000001200_df_long <- F_somatic_ENSPREG00000001200_df_long [order(F_somatic_ENSPREG00000001200_df_long$count, decreasing = F), ]

F_somatic_ENSPREG00000005616_df <- somatic_allele_counts %>% filter(geneid == "ENSPREG00000005616" & sex == "F")
F_somatic_ENSPREG00000005616_df_long <- F_somatic_ENSPREG00000005616_df %>%
  pivot_longer(cols = c(ref_count, alt_count), names_to = "count_type", values_to = "count") %>%
  mutate(count_type = case_when(count_type == "ref_count" ~ "Reference", count_type == "alt_count" ~ "Alternate"))
F_somatic_ENSPREG00000005616_df_long <- F_somatic_ENSPREG00000005616_df_long [order(F_somatic_ENSPREG00000005616_df_long$count, decreasing = F), ]

F_somatic_ENSPREG00000019478_df <- somatic_allele_counts %>% filter(geneid == "ENSPREG00000019478" & sex == "F")
F_somatic_ENSPREG00000019478_df_long <- F_somatic_ENSPREG00000019478_df %>%
  pivot_longer(cols = c(ref_count, alt_count), names_to = "count_type", values_to = "count") %>%
  mutate(count_type = case_when(count_type == "ref_count" ~ "Reference", count_type == "alt_count" ~ "Alternate"))
F_somatic_ENSPREG00000019478_df_long <- F_somatic_ENSPREG00000019478_df_long [order(F_somatic_ENSPREG00000019478_df_long$count, decreasing = F), ]

F_somatic_unc5da_df <- somatic_allele_counts %>% filter(geneid == "unc5da" & sex == "F")
F_somatic_unc5da_df_long <- F_somatic_unc5da_df %>%
  pivot_longer(cols = c(ref_count, alt_count), names_to = "count_type", values_to = "count") %>%
  mutate(count_type = case_when(count_type == "ref_count" ~ "Reference", count_type == "alt_count" ~ "Alternate"))
F_somatic_unc5da_df_long <- F_somatic_unc5da_df_long [order(F_somatic_unc5da_df_long$count, decreasing = F), ]

wilcox.test(F_somatic_CACNA1B_df_long$count ~ F_somatic_CACNA1B_df_long$count_type, p.adjust.methods = "bonferroni") # 0.09376
wilcox.test(F_somatic_ENSPREG00000001200_df_long$count ~ F_somatic_ENSPREG00000001200_df_long$count_type, p.adjust.methods = "bonferroni") # 0.003104
wilcox.test(F_somatic_ENSPREG00000005616_df_long$count ~ F_somatic_ENSPREG00000005616_df_long$count_type, p.adjust.methods = "bonferroni") # 0.4
wilcox.test(F_somatic_ENSPREG00000019478_df_long$count ~ F_somatic_ENSPREG00000019478_df_long$count_type, p.adjust.methods = "bonferroni") # 0.4111
wilcox.test(F_somatic_unc5da_df_long$count ~ F_somatic_unc5da_df_long$count_type, p.adjust.methods = "bonferroni") # 0.01421


#### Boxplot ####

M_somatic_genes_combo <- rbind(M_somatic_CACNA1B_df_long, M_somatic_ENSPREG00000001200_df_long, M_somatic_ENSPREG00000005616_df_long,
                               M_somatic_ENSPREG00000019478_df_long, M_somatic_unc5da_df_long)
M_somatic_genes_combo <- M_somatic_genes_combo %>% mutate(geneid = case_when(
  geneid == "unc5da" ~ "UNC5DA", 
  geneid == "ENSPREG00000001200" ~ "SPP1", 
  geneid == "ENSPREG00000019478" ~ "PUM3", 
  TRUE ~ geneid))

plot_somatic_boxplot <- ggplot(M_somatic_genes_combo, aes(x = geneid, y = count, fill = count_type, col = count_type)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 4, lwd = 1.5) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 65)) +
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 4, alpha = 0.5) +
  scale_fill_manual(values = c("seagreen", "slategray3")) + scale_colour_manual(values = c("darkgreen", "slategray4")) +
  ylab("Allele Count") + xlab("Gene") + labs(fill = "Allele Type", col = "Allele Type") +
  theme_classic() + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), 
                          axis.text.x = element_text(size = 16, angle = 25, vjust = 0.95, hjust = 0.9, face = "bold"), 
                          axis.ticks.x = element_blank(), legend.position = "top",
                          axis.text.y = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5),
                          plot.margin = margin(t = 20, r = 50, b = 50, l = 20))

plot_somatic_boxplot


F_somatic_genes_combo <- rbind(F_somatic_CACNA1B_df_long, F_somatic_ENSPREG00000001200_df_long, F_somatic_ENSPREG00000005616_df_long,
                               F_somatic_ENSPREG00000019478_df_long, F_somatic_unc5da_df_long)
F_somatic_genes_combo <- F_somatic_genes_combo %>% mutate(geneid = case_when(
  geneid == "unc5da" ~ "UNC5DA", 
  geneid == "ENSPREG00000001200" ~ "SPP1", 
  geneid == "ENSPREG00000019478" ~ "PUM3", 
  TRUE ~ geneid))

F_plot_somatic_boxplot <- ggplot(F_somatic_genes_combo, aes(x = geneid, y = count, fill = count_type, col = count_type)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 4, lwd = 1.5) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 65)) +
  scale_fill_manual(values = c("pink", "slategray3")) + scale_colour_manual(values = c("pink3", "slategray4")) +
  ylab("Allele Count") + xlab("Gene") + labs(fill = "Allele Type", col = "Allele Type") + ggtitle("") +
  theme_classic() + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
                          axis.text.x = element_text(size = 16, angle = 25, vjust = 0.95, hjust = 0.9, face = "bold"), 
                          axis.ticks.x = element_blank(),legend.position = "top",
                          axis.text.y = element_text(size = 18), plot.title = element_text(size = 25, hjust = 0.5),
                          plot.margin = margin(t = 10, r = 10, b = 10, l = 10))

F_plot_somatic_boxplot

#### Now gonadal tissues ####

M_gonad_CACNA1B_df <- gonad_allele_counts %>% filter(geneid == "CACNA1B" & sex == "M")
M_gonad_CACNA1B_df$total_count <- M_gonad_CACNA1B_df$ref_count + M_gonad_CACNA1B_df$alt_count
M_gonad_CACNA1B_df_long <- M_gonad_CACNA1B_df %>%
  pivot_longer(cols = c(ref_count, alt_count), names_to = "count_type", values_to = "count") %>%
  mutate(count_type = case_when(count_type == "ref_count" ~ "Reference", count_type == "alt_count" ~ "Alternate"))
M_gonad_CACNA1B_df_long <- M_gonad_CACNA1B_df_long [order(M_gonad_CACNA1B_df_long$count, decreasing = F), ]

M_gonad_ENSPREG00000001200_df <- gonad_allele_counts %>% filter(geneid == "ENSPREG00000001200" & sex == "M")
M_gonad_ENSPREG00000001200_df$total_count <- M_gonad_ENSPREG00000001200_df$ref_count + M_gonad_ENSPREG00000001200_df$alt_count
M_gonad_ENSPREG00000001200_df_long <- M_gonad_ENSPREG00000001200_df %>%
  pivot_longer(cols = c(ref_count, alt_count), names_to = "count_type", values_to = "count") %>%
  mutate(count_type = case_when(count_type == "ref_count" ~ "Reference", count_type == "alt_count" ~ "Alternate"))
M_gonad_ENSPREG00000001200_df_long <- M_gonad_ENSPREG00000001200_df_long [order(M_gonad_ENSPREG00000001200_df_long$count, decreasing = F), ]

M_gonad_ENSPREG00000005616_df <- gonad_allele_counts %>% filter(geneid == "ENSPREG00000005616" & sex == "M")
M_gonad_ENSPREG00000005616_df$total_count <- M_gonad_ENSPREG00000005616_df$ref_count + M_gonad_ENSPREG00000005616_df$alt_count
M_gonad_ENSPREG00000005616_df_long <- M_gonad_ENSPREG00000005616_df %>%
  pivot_longer(cols = c(ref_count, alt_count), names_to = "count_type", values_to = "count") %>%
  mutate(count_type = case_when(count_type == "ref_count" ~ "Reference", count_type == "alt_count" ~ "Alternate"))
M_gonad_ENSPREG00000005616_df_long <- M_gonad_ENSPREG00000005616_df_long [order(M_gonad_ENSPREG00000005616_df_long$count, decreasing = F), ]

M_gonad_ENSPREG00000019478_df <- gonad_allele_counts %>% filter(geneid == "ENSPREG00000019478" & sex == "M")
M_gonad_ENSPREG00000019478_df$total_count <- M_gonad_ENSPREG00000019478_df$ref_count + M_gonad_ENSPREG00000019478_df$alt_count
M_gonad_ENSPREG00000019478_df_long <- M_gonad_ENSPREG00000019478_df %>%
  pivot_longer(cols = c(ref_count, alt_count), names_to = "count_type", values_to = "count") %>%
  mutate(count_type = case_when(count_type == "ref_count" ~ "Reference", count_type == "alt_count" ~ "Alternate"))
M_gonad_ENSPREG00000019478_df_long <- M_gonad_ENSPREG00000019478_df_long [order(M_gonad_ENSPREG00000019478_df_long$count, decreasing = F), ]

M_gonad_unc5da_df <- gonad_allele_counts %>% filter(geneid == "unc5da" & sex == "M")
M_gonad_unc5da_df$total_count <- M_gonad_unc5da_df$ref_count + M_gonad_unc5da_df$alt_count
M_gonad_unc5da_df_long <- M_gonad_unc5da_df %>%
  pivot_longer(cols = c(ref_count, alt_count), names_to = "count_type", values_to = "count") %>%
  mutate(count_type = case_when(count_type == "ref_count" ~ "Reference", count_type == "alt_count" ~ "Alternate"))
M_gonad_unc5da_df_long <- M_gonad_unc5da_df_long [order(M_gonad_unc5da_df_long$count, decreasing = F), ]

wilcox.test(M_gonad_CACNA1B_df_long$count ~ M_gonad_CACNA1B_df_long$count_type, p.adjust.methods = "bonferroni") # 0.2477
wilcox.test(M_gonad_ENSPREG00000001200_df_long$count ~ M_gonad_ENSPREG00000001200_df_long$count_type, p.adjust.methods = "bonferroni") # 0.1882
wilcox.test(M_gonad_unc5da_df_long$count ~ M_gonad_unc5da_df_long$count_type, p.adjust.methods = "bonferroni") # 0.1617

#### For Female ####

F_gonad_CACNA1B_df <- gonad_allele_counts %>% filter(geneid == "CACNA1B" & sex == "F")
F_gonad_CACNA1B_df_long <- F_gonad_CACNA1B_df %>%
  pivot_longer(cols = c(ref_count, alt_count), names_to = "count_type", values_to = "count") %>%
  mutate(count_type = case_when(count_type == "ref_count" ~ "Reference", count_type == "alt_count" ~ "Alternate"))
F_gonad_CACNA1B_df_long <- F_gonad_CACNA1B_df_long [order(F_gonad_CACNA1B_df_long$count, decreasing = F), ]

F_gonad_ENSPREG00000001200_df <- gonad_allele_counts %>% filter(geneid == "ENSPREG00000001200" & sex == "F")
F_gonad_ENSPREG00000001200_df_long <- F_gonad_ENSPREG00000001200_df %>%
  pivot_longer(cols = c(ref_count, alt_count), names_to = "count_type", values_to = "count") %>%
  mutate(count_type = case_when(count_type == "ref_count" ~ "Reference", count_type == "alt_count" ~ "Alternate"))
F_gonad_ENSPREG00000001200_df_long <- F_gonad_ENSPREG00000001200_df_long [order(F_gonad_ENSPREG00000001200_df_long$count, decreasing = F), ]

F_gonad_ENSPREG00000005616_df <- gonad_allele_counts %>% filter(geneid == "ENSPREG00000005616" & sex == "F")
F_gonad_ENSPREG00000005616_df_long <- F_gonad_ENSPREG00000005616_df %>%
  pivot_longer(cols = c(ref_count, alt_count), names_to = "count_type", values_to = "count") %>%
  mutate(count_type = case_when(count_type == "ref_count" ~ "Reference", count_type == "alt_count" ~ "Alternate"))
F_gonad_ENSPREG00000005616_df_long <- F_gonad_ENSPREG00000005616_df_long [order(F_gonad_ENSPREG00000005616_df_long$count, decreasing = F), ]

F_gonad_ENSPREG00000019478_df <- gonad_allele_counts %>% filter(geneid == "ENSPREG00000019478" & sex == "F")
F_gonad_ENSPREG00000019478_df_long <- F_gonad_ENSPREG00000019478_df %>%
  pivot_longer(cols = c(ref_count, alt_count), names_to = "count_type", values_to = "count") %>%
  mutate(count_type = case_when(count_type == "ref_count" ~ "Reference", count_type == "alt_count" ~ "Alternate"))
F_gonad_ENSPREG00000019478_df_long <- F_gonad_ENSPREG00000019478_df_long [order(F_gonad_ENSPREG00000019478_df_long$count, decreasing = F), ]

F_gonad_unc5da_df <- gonad_allele_counts %>% filter(geneid == "unc5da" & sex == "F")
F_gonad_unc5da_df_long <- F_gonad_unc5da_df %>%
  pivot_longer(cols = c(ref_count, alt_count), names_to = "count_type", values_to = "count") %>%
  mutate(count_type = case_when(count_type == "ref_count" ~ "Reference", count_type == "alt_count" ~ "Alternate"))
F_gonad_unc5da_df_long <- F_gonad_unc5da_df_long [order(F_gonad_unc5da_df_long$count, decreasing = F), ]

wilcox.test(F_gonad_CACNA1B_df_long$count ~ F_gonad_CACNA1B_df_long$count_type, p.adjust.methods = "bonferroni") # 0.01442
wilcox.test(F_gonad_unc5da_df_long$count ~ F_gonad_unc5da_df_long$count_type, p.adjust.methods = "bonferroni") # 0.1547


#### Gonad Barplot ####

M_gonad_genes_combo <- rbind(M_gonad_CACNA1B_df_long, M_gonad_ENSPREG00000001200_df_long, M_gonad_unc5da_df_long)
M_gonad_genes_combo <- M_gonad_genes_combo %>% mutate(geneid = case_when(
  geneid == "unc5da" ~ "UNC5DA", 
  geneid == "ENSPREG00000001200" ~ "SPP1",
  TRUE ~ geneid))

plot_gonad_boxplot <- ggplot(M_gonad_genes_combo, aes(x = geneid, y = count, fill = count_type, col = count_type)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 4, lwd = 1.5) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 65)) +
  scale_fill_manual(values = c("seagreen", "slategray3")) + scale_colour_manual(values = c("darkgreen", "slategray4")) +
  ylab("Allele Count") + xlab("Gene") + labs(fill = "Allele Type", col = "Allele Type") +
  theme_classic() + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), 
                          axis.text.x = element_text(size = 16, angle = 25, vjust = 0.95, hjust = 0.9, face = "bold"), 
                          axis.ticks.x = element_blank(), legend.position = "top",
                          axis.text.y = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5),
                          plot.margin = margin(t = 20, r = 50, b = 50, l = 20))

plot_gonad_boxplot

F_gonad_genes_combo <- rbind(F_gonad_CACNA1B_df_long, F_gonad_unc5da_df_long)
F_gonad_genes_combo <- F_gonad_genes_combo %>% mutate(geneid = case_when(geneid == "unc5da" ~ "UNC5DA", TRUE ~ geneid))

F_plot_gonad_boxplot <- ggplot(F_gonad_genes_combo, aes(x = geneid, y = count, fill = count_type, col = count_type)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 4, lwd = 1.5) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 65)) +
  scale_fill_manual(values = c("pink", "slategray3")) + scale_colour_manual(values = c("pink3", "slategray4")) +
  ylab("Allele Count") + xlab("Gene") + labs(fill = "Allele Type", col = "Allele Type") + ggtitle("") +
  theme_classic() + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), 
                          axis.text.x = element_text(size = 16, angle = 25, vjust = 0.95, hjust = 0.9, face = "bold"), 
                          axis.ticks.x = element_blank(), legend.position = "top",
                          axis.text.y = element_text(size = 18), plot.title = element_text(size = 25, hjust = 0.5),
                          plot.margin = margin(t = 20, r = 50, b = 50, l = 20))

F_plot_gonad_boxplot


### Final collection of plots:

(plot_gonad_A | plot_gonad_SC | plot_gonad_boxplot) / (plot_somatic_A | plot_somatic_SC | plot_somatic_boxplot) +
  plot_layout(guides = "collect")

(F_plot_gonad_boxplot | F_plot_somatic_boxplot) + plot_layout(guides = "collect")

### Version 2 of plots:

plot_somatic_boxplot_v2 <- ggplot(M_somatic_genes_combo, aes(x = geneid, y = count, fill = count_type, col = count_type)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 4, lwd = 1.5) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 65)) +
  scale_fill_manual(values = c("seagreen", "slategray3")) + scale_colour_manual(values = c("darkgreen", "slategray4")) +
  ylab("Allele Count") + xlab("") + labs(fill = "Allele Type", col = "Allele Type") +
  theme_classic() + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), 
                          axis.text.x = element_text(size = 16, angle = 25, vjust = 0.95, hjust = 0.9, face = "bold"), 
                          axis.ticks.x = element_blank(), legend.position = "top",
                          legend.title = element_text(size = 18), legend.text = element_text(size = 18),
                          axis.text.y = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5),
                          plot.margin = margin(b = 0))

plot_gonad_boxplot_v2 <- ggplot(M_gonad_genes_combo, aes(x = geneid, y = count, fill = count_type, col = count_type)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 4, lwd = 1.5) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 65)) +
  scale_fill_manual(values = c("seagreen", "slategray3")) + scale_colour_manual(values = c("darkgreen", "slategray4")) +
  ylab("") + xlab("") + labs(fill = "Allele Type", col = "Allele Type") +
  theme_classic() + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), 
                          axis.text.x = element_text(size = 16, angle = 25, vjust = 0.95, hjust = 0.9, face = "bold"), 
                          axis.ticks.x = element_blank(), legend.position = "top",
                          legend.title = element_text(size = 18), legend.text = element_text(size = 18),
                          axis.text.y = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5),
                          plot.margin = margin(b = 0))

F_plot_somatic_boxplot_v2 <- ggplot(F_somatic_genes_combo, aes(x = geneid, y = count, fill = count_type, col = count_type)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 4, lwd = 1.5) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 65)) +
  scale_fill_manual(values = c("pink", "slategray3")) + scale_colour_manual(values = c("pink3", "slategray4")) +
  ylab("Allele Count") + xlab("Gene") + labs(fill = "Allele Type", col = "Allele Type") + ggtitle("") +
  theme_classic() + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
                          axis.text.x = element_text(size = 16, angle = 25, vjust = 0.95, hjust = 0.9, face = "bold"), 
                          axis.ticks.x = element_blank(),legend.position = "top",
                          legend.title = element_text(size = 18), legend.text = element_text(size = 18),
                          axis.text.y = element_text(size = 18), plot.title = element_text(size = 25, hjust = 0.5),
                          plot.margin = margin(t = 0))

F_plot_gonad_boxplot_v2 <- ggplot(F_gonad_genes_combo, aes(x = geneid, y = count, fill = count_type, col = count_type)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 4, lwd = 1.5) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 65)) +
  scale_fill_manual(values = c("pink", "slategray3")) + scale_colour_manual(values = c("pink3", "slategray4")) +
  ylab("") + xlab("Gene") + labs(fill = "Allele Type", col = "Allele Type") + ggtitle("") +
  theme_classic() + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), 
                          axis.text.x = element_text(size = 16, angle = 25, vjust = 0.95, hjust = 0.9, face = "bold"), 
                          axis.ticks.x = element_blank(), legend.position = "top",
                          legend.title = element_text(size = 18), legend.text = element_text(size = 18),
                          axis.text.y = element_text(size = 18), plot.title = element_text(size = 25, hjust = 0.5),
                          plot.margin = margin(t = 0))


(plot_somatic_boxplot_v2 | plot_gonad_boxplot_v2)/(F_plot_somatic_boxplot_v2 | F_plot_gonad_boxplot_v2) 

clean_plot_somatic_male_v2 <- ggplot(clean_male_somatic_allele_counts, aes(x = mar)) + geom_density(aes(fill = chrom_type, col = chrom_type), alpha = 0.6, adjust = 1.5, linewidth = 1.2) + 
  scale_fill_manual(values = c("#857C7B", "mediumpurple3")) + labs(fill = "SNP Type", col = "SNP Type") +
  scale_colour_manual(values = c("#857C7B", "mediumpurple3")) + 
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +  scale_y_continuous(expand = c(0,0), lim = c(0, 7)) +
  ylab("Density") + xlab("") + ggtitle("") +
  geom_vline(aes(xintercept=0.6071429),color="grey40",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.617),color="darkorchid3",linetype="dashed",size=0.8) +
  theme_classic() + theme(
    axis.title.x = element_text(size = 25),  axis.title.y = element_text(size = 25), legend.position = "top",
    axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), 
    legend.title = element_text(size = 18), legend.text = element_text(size = 18),
    plot.title = element_text(size = 20, hjust = 0.5),
    plot.margin = margin(r = 20)
  )

clean_plot_gonad_male_v2 <- ggplot(clean_male_gonad_allele_counts, aes(x = mar)) + geom_density(aes(fill = chrom_type, col = chrom_type), alpha = 0.5, adjust = 1.5, linewidth = 1.2) + 
  scale_fill_manual(values = c("#857C7B", "mediumpurple3")) + labs(fill = "SNP Type", col = "SNP Type") +
  scale_colour_manual(values = c("#857C7B", "mediumpurple3")) + 
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +  scale_y_continuous(expand = c(0,0), lim = c(0, 7)) +
  ylab("") + xlab("") + ggtitle("") +
  geom_vline(aes(xintercept=0.5862069),color="grey40",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5777778),color="darkorchid3",linetype="dashed",size=0.8) +
  theme_classic() + theme(
    axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25), legend.position = "top",
    axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
    legend.title = element_text(size = 18), legend.text = element_text(size = 18),
    plot.title = element_text(size = 20, hjust = 0.5),
    plot.margin = margin(r = 20)
  )

clean_plot_somatic_female_v2 <- ggplot(clean_female_somatic_allele_counts, aes(x = mar)) + geom_density(aes(fill = chrom_type, col = chrom_type), alpha = 0.5, adjust = 1.5, linewidth = 1.2) + 
  scale_fill_manual(values = c("#857C7B", "mediumpurple3")) + labs(fill = "SNP Type", col = "SNP Type") +
  scale_colour_manual(values = c("#857C7B", "mediumpurple3")) + 
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +  scale_y_continuous(expand = c(0,0), lim = c(0, 7)) +
  ylab("Density") + xlab("Major Allele Ratio") + ggtitle("") +
  geom_vline(aes(xintercept=0.5862069),color="grey40",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5714286),color="darkorchid3",linetype="dashed",size=0.8) +
  theme_classic() + theme(
    axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25), legend.position = "none",
    axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
    plot.title = element_text(size = 20, hjust = 0.5),
    plot.margin = margin(r = 20)
  )

clean_plot_gonad_female_v2 <- ggplot(clean_female_gonad_allele_counts, aes(x = mar)) + geom_density(aes(fill = chrom_type, col = chrom_type), alpha = 0.5, adjust = 1.5, linewidth = 1.2) + 
  scale_fill_manual(values = c("#857C7B", "mediumpurple3")) + labs(fill = "SNP Type", col = "SNP Type") +
  scale_colour_manual(values = c("#857C7B", "mediumpurple3")) + 
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +  scale_y_continuous(expand = c(0,0), lim = c(0, 7)) +
  ylab("") + xlab("Major Allele Ratio") + ggtitle("") +
  geom_vline(aes(xintercept=0.5957447),color="grey40",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6060606),color="darkorchid3",linetype="dashed",size=0.8) +
  theme_classic() + theme(
    axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25), legend.position = "none",
    axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
    plot.title = element_text(size = 20, hjust = 0.5),
    plot.margin = margin(r = 20)
  )


(clean_plot_somatic_male_v2 | clean_plot_gonad_male_v2) / (clean_plot_somatic_female_v2 | clean_plot_gonad_female_v2)


