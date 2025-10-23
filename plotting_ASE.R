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

gonad_minfilt <- read.table(file = "gonad_filtered10.txt", header = T)
skin_minfilt <- read.table(file = "skin_filtered10.txt", header = T)
liver_minfilt <- read.table(file = "liver_filtered10.txt", header = T)
heart_minfilt <- read.table(file = "heart_filtered10.txt", header = T)

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

somatic_all %>% pull(alt_count) %>% quantile(probs = 0.8) 
somatic_all %>% pull(ref_count) %>% quantile(probs = 0.8)
gonad_minfilt %>% pull(alt_count) %>% quantile(probs = 0.8)
gonad_minfilt %>% pull(ref_count) %>% quantile(probs = 0.8)

somatic_allele_counts <- somatic_all %>%  filter(ref_count <= 49 & alt_count <= 104)
gonad_allele_counts <- gonad_all %>%  filter(ref_count <= 41 & alt_count <= 78)

#################################################################
########### Organize tissues, then plot and stats ###############
#################################################################

gonad_ase <- gonad_allele_counts
gonad_ase$total_counts <- gonad_ase$ref_count + gonad_ase$alt_count
gonad_ase$mar <- pmax(gonad_ase$ref_count, gonad_ase$alt_count)/(gonad_ase$total_counts)
gonad_ase <- na.omit(gonad_ase)

gonad_samp_filt <- gonad_ase %>% group_by(geneid) %>%
  filter(sum(sex == "F") >= 3, sum(sex == "M") >= 3) %>% ungroup()

somatic_ase <- somatic_allele_counts
somatic_ase$total_counts <- somatic_ase$ref_count + somatic_ase$alt_count
somatic_ase$mar <- pmax(somatic_ase$ref_count, somatic_ase$alt_count)/(somatic_ase$total_counts)
somatic_ase <- na.omit(somatic_ase)

somatic_samp_filt <- somatic_ase %>% group_by(geneid) %>%
  filter(sum(sex == "F") >= 3, sum(sex == "M") >= 3) %>% ungroup()

#############################################################################
################### Plot all impacted gene allele counts: ####################
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

wilcox.test(M_somatic_CACNA1B_df_long$count ~ M_somatic_CACNA1B_df_long$count_type, p.adjust.methods = "bonferroni")
wilcox.test(M_somatic_ENSPREG00000001200_df_long$count ~ M_somatic_ENSPREG00000001200_df_long$count_type, p.adjust.methods = "bonferroni")
wilcox.test(M_somatic_ENSPREG00000005616_df_long$count ~ M_somatic_ENSPREG00000005616_df_long$count_type, p.adjust.methods = "bonferroni") 
wilcox.test(M_somatic_ENSPREG00000019478_df_long$count ~ M_somatic_ENSPREG00000019478_df_long$count_type, p.adjust.methods = "bonferroni") 
wilcox.test(M_somatic_unc5da_df_long$count ~ M_somatic_unc5da_df_long$count_type, p.adjust.methods = "bonferroni")

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

wilcox.test(F_somatic_CACNA1B_df_long$count ~ F_somatic_CACNA1B_df_long$count_type, p.adjust.methods = "bonferroni")
wilcox.test(F_somatic_ENSPREG00000001200_df_long$count ~ F_somatic_ENSPREG00000001200_df_long$count_type, p.adjust.methods = "bonferroni") 
wilcox.test(F_somatic_ENSPREG00000005616_df_long$count ~ F_somatic_ENSPREG00000005616_df_long$count_type, p.adjust.methods = "bonferroni")
wilcox.test(F_somatic_ENSPREG00000019478_df_long$count ~ F_somatic_ENSPREG00000019478_df_long$count_type, p.adjust.methods = "bonferroni")
wilcox.test(F_somatic_unc5da_df_long$count ~ F_somatic_unc5da_df_long$count_type, p.adjust.methods = "bonferroni") 

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

wilcox.test(M_gonad_CACNA1B_df_long$count ~ M_gonad_CACNA1B_df_long$count_type, p.adjust.methods = "bonferroni")
wilcox.test(M_gonad_ENSPREG00000001200_df_long$count ~ M_gonad_ENSPREG00000001200_df_long$count_type, p.adjust.methods = "bonferroni") 
wilcox.test(M_gonad_unc5da_df_long$count ~ M_gonad_unc5da_df_long$count_type, p.adjust.methods = "bonferroni") 

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

wilcox.test(F_gonad_CACNA1B_df_long$count ~ F_gonad_CACNA1B_df_long$count_type, p.adjust.methods = "bonferroni")
wilcox.test(F_gonad_unc5da_df_long$count ~ F_gonad_unc5da_df_long$count_type, p.adjust.methods = "bonferroni")


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


