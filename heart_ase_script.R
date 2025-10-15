###################################################
############ Using Iulia's scRNA seq ##############
############ To identify ASE in genes #############
###################################################

rm(list=ls())
ls() 

setwd("//files.zoology.ubc.ca/ljmfong/flex")
getwd()

#Ensure you're on the latest R version

###################################################
########## 1. Load in required packages ###########
###################################################

library(cowplot)
library(DESeq2)
library(dplyr)
library(ggbreak)
library(ggplot2)
library(magrittr)
library(Matrix)
library(Matrix.utils)
library(purrr)
library(RCurl)
library(scales)
library(scuttle)
library(Seurat)
library(SingleCellExperiment)
library(stringr)
library(tibble)
library(tidyverse)
library(patchwork)

###################################################
########### 2. Load in Iulia's matrices ###########
###################################################

load("iulia_scRNA/heart_umap_3_rename_clusters_dbrmv.RData")

counts_heart <- heart_umap_3_rename_clusters_dbrmv@assays$RNA@counts
metadata_heart <- heart_umap_3_rename_clusters_dbrmv@meta.data
metadata_heart$CellType <- factor(heart_umap_3_rename_clusters_dbrmv@active.ident)

#### 3. Get differential gene expression for separate cell types ####

#allelecount_Lib5 <- read.table("scrna_scAlleleCount/complete_matrice_for_R/Lib5_counts.txt", sep = " ", header = T)
#allelecount_Lib16 <- read.table("scrna_scAlleleCount/complete_matrice_for_R/Lib16_counts.txt", sep = " ", header = T)
#allelecount_Lib27 <- read.table("scrna_scAlleleCount/complete_matrice_for_R/Lib27_counts.txt", sep = " ", header = T)
#allelecount_Lib6 <- read.table("scrna_scAlleleCount/complete_matrice_for_R/Lib6_counts.txt", sep = " ", header = T)
#allelecount_Lib7 <- read.table("scrna_scAlleleCount/complete_matrice_for_R/Lib7_counts.txt", sep = " ", header = T)
#allelecount_Lib21 <- read.table("scrna_scAlleleCount/complete_matrice_for_R/Lib21_counts.txt", sep = " ", header = T)

allelecount_Lib5 <- read.table("cellSNP_scAlleleCount/complete_matrice_for_R/Lib5_counts.txt", sep = " ", header = T)
allelecount_Lib16 <- read.table("cellSNP_scAlleleCount/complete_matrice_for_R/Lib16_counts.txt", sep = " ", header = T)
allelecount_Lib27 <- read.table("cellSNP_scAlleleCount/complete_matrice_for_R/Lib27_counts.txt", sep = " ", header = T)
allelecount_Lib6 <- read.table("cellSNP_scAlleleCount/complete_matrice_for_R/Lib6_counts.txt", sep = " ", header = T)
allelecount_Lib7 <- read.table("cellSNP_scAlleleCount/complete_matrice_for_R/Lib7_counts.txt", sep = " ", header = T)
allelecount_Lib21 <- read.table("cellSNP_scAlleleCount/complete_matrice_for_R/Lib21_counts.txt", sep = " ", header = T)

merged_heart_AC <- rbind(allelecount_Lib5, allelecount_Lib16, allelecount_Lib27, 
                        allelecount_Lib6, allelecount_Lib7, allelecount_Lib21)
colnames(merged_heart_AC)[3] <- "cells"

#### Attach gene info and find matching gene_name for positional information:

gene_info <- read.table("new_Ret_SNPs/gene_boundary_w_gene_names.bed", header = T, sep = "\t")

attached_gene_AC <- merged_heart_AC %>%
  mutate(CHROM = str_extract(j, "^[^:]+"),
         POS = as.numeric(str_extract(j, "\\d+(?=-)"))
  ) # Extracts CHROM and POS from the 'j' column

#Find the matching gene_name for the positional info:

attached_gene_AC <- attached_gene_AC %>%
  rowwise() %>%
  mutate(
    geneid = {
      match_row <- gene_info %>%
        filter(CHROM == CHROM & START <= POS & END >= POS)
      if (nrow(match_row) > 0) match_row$GENE_NAME[1] else NA
    }
  ) %>%
  ungroup()

cleaned_df <- select(attached_gene_AC, "cells", "ref", "alt", "CHROM", "POS", "geneid")
#write.table(cleaned_df, "scASE_results/heart_cleaned_df_gene_attached.txt", quote = F, row.names = F)
write.table(cleaned_df, "cellSNP_scAlleleCount/R_results/heart_cleaned_df_gene_attached.txt", quote = F, row.names = F)
#If loading in:
#cleaned_df <- read.table(file = "scASE_results/heart_cleaned_df_gene_attached.txt", header = T)
cleaned_df <- read.table(file = "cellSNP_scAlleleCount/R_results/heart_cleaned_df_gene_attached.txt", header = T)

######## Attach to your SingleCellExperiment ########

#OG sce_heart command:
sce_heart <- SingleCellExperiment(assays=list(counts=counts_heart), colData=metadata_heart)

# Determine the number of clusters (nk) and the cluster names (kids), and the number of samples (ns) and sample names (sids)
colData(sce_heart)$sample_id <- as.factor(colData(sce_heart)$sample)
colData(sce_heart)$sex_id <- as.factor(colData(sce_heart)$sex)
colData(sce_heart)$cluster_id <- as.factor(colData(sce_heart)$CellType)
kids_heart <- purrr::set_names(levels(sce_heart$cluster_id))
nk_heart <- length(kids_heart)
sids_heart <- purrr::set_names(levels(sce_heart$sample_id))
ns_heart <- length(sids_heart)

# Turn named vector into a numeric vector
n_cells_heart <- as.numeric(table(sce_heart$sample_id))

# Reorder samples (rows) of the metadata to match the order of the sample names
m_heart <- match(sids_heart, sce_heart$sample_id)

# Create sample level metadata by combining the reordered metadata with the number of cells
ei_heart <- data.frame(colData(sce_heart)[m_heart, ], n_cells_heart, row.names=NULL) %>% select(-"cluster_id")

# Remove lowly expressed genes
sce_heart <- sce_heart[rowSums(counts(sce_heart)>1)>=10,]	

# Aggregate counts per sample_id and cluster_id
groups_heart <- colData(sce_heart)[, c("cluster_id", "sample_id")]
pb_heart <- aggregate.Matrix(t(counts(sce_heart)), groupings=groups_heart, fun="sum")


##############

summed_ac_counts <- cleaned_df %>%
  group_by(cells, geneid) %>%
  summarise(ref_count = sum(ref, na.rm = T), alt_count = sum(alt, na.rm = T), .groups = "drop")

any(duplicated(summed_ac_counts$cells)) #This should be TRUE

summed_ac_counts <- summed_ac_counts %>%
  mutate(cells = str_replace(cells, "^Lib6", "F1_Lib6"),
         cells = str_replace(cells, "^Lib7", "F2_Lib7"),
         cells = str_replace(cells, "^Lib21", "F3_Lib21"),
         cells = str_replace(cells, "^Lib5", "M1_Lib5"),
         cells = str_replace(cells, "^Lib16", "M2_Lib16"),
         cells = str_replace(cells, "^Lib27", "M3_Lib27"))

### Cluster your summed_ac_counts into cell-types
#All the info for the match is in groups_hearts

groups_heart_name <- data.frame(groups_heart) 
groups_heart_name <- groups_heart_name %>% rownames_to_column(var = "cells")

match(groups_heart_name$cells, summed_ac_counts$cells)

group_name_ac_heart <- summed_ac_counts %>%
  left_join(groups_heart_name, by = "cells")

cleaned_group_name_ac_heart <- na.omit(group_name_ac_heart)

################ Turn into matrix ####################
############## save as dgCMatrix #####################

# Build ONE big sparse matrix (cells x genes) for ref_count & alt_count
ref_matrix_all <- cleaned_group_name_ac_heart %>%
  select(cells, geneid, ref_count) %>%
  group_by(cells, geneid) %>%
  summarise(count = sum(ref_count), .groups = "drop") %>%
  pivot_wider(names_from = geneid, values_from = count, values_fill = 0) %>%
  column_to_rownames("cells") %>%
  as.matrix() %>%
  Matrix(sparse = TRUE)


alt_matrix_all <- cleaned_group_name_ac_heart %>%
  select(cells, geneid, alt_count) %>%
  group_by(cells, geneid) %>%
  summarise(count = sum(alt_count), .groups = "drop") %>%
  pivot_wider(names_from = geneid, values_from = count, values_fill = 0) %>%
  column_to_rownames("cells") %>%
  as.matrix() %>%
  Matrix(sparse = TRUE)

cell_metadata <- cleaned_group_name_ac_heart %>%
  select(cells, cluster_id, sample_id) %>%
  distinct()

# Ensure the rownames of matrix match metadata$cell
ref_matrix_all <- ref_matrix_all[rownames(ref_matrix_all) %in% cell_metadata$cells, ]
cell_metadata <- cell_metadata[match(rownames(ref_matrix_all), cell_metadata$cells), ]

# Combine group labels
cell_metadata$group <- paste(cell_metadata$cluster_id, cell_metadata$sample_id, sep = "_")

# Aggregate by group
aggregated_ref <- aggregate.Matrix(ref_matrix_all, groupings = cell_metadata$group, fun = "sum")
aggregated_alt <- aggregate.Matrix(alt_matrix_all, groupings = cell_metadata$group, fun = "sum")

###################################################
# Match matrices style to Iulia's

ref_splitf_heart <- sapply(stringr::str_split(rownames(aggregated_ref), pattern="_", n=2), '[', 1) # Check that each cluster is present in all samples
alt_splitf_heart <- sapply(stringr::str_split(rownames(aggregated_alt), pattern="_", n=2), '[', 1) # Check that each cluster is present in all samples

ref_pb_heart <- split.data.frame(aggregated_ref, factor(ref_splitf_heart)) %>% lapply(function(u) set_colnames(t(u), stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))
alt_pb_heart <- split.data.frame(aggregated_alt, factor(alt_splitf_heart)) %>% lapply(function(u) set_colnames(t(u), stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

##############################################################################
########### Find MAR on all celltypes before splitting them up: ##############
##############################################################################

working_ac_heart <- summed_ac_counts
working_ac_heart <- working_ac_heart %>%
  filter(ref_count >= 10 & alt_count >= 10)
working_ac_heart$log10_ref <- log10(working_ac_heart$ref_count)
working_ac_heart$log10_alt <- log10(working_ac_heart$alt_count)

working_ac_heart <- working_ac_heart %>%
  mutate(sample_id = case_when( 
           grepl("^M1", cells) ~ "M1",
           grepl("^M2", cells) ~ "M2",
           grepl("^M3", cells) ~ "M3",
           grepl("^F1", cells) ~ "F1",
           grepl("^F2", cells) ~ "F2",
           grepl("^F3", cells) ~ "F3"))
#write.table(working_ac_heart, "scASE_results/heart_filtered10.txt", quote = F, row.names = F)
write.table(working_ac_heart, "cellSNP_scAlleleCount/R_results/heart_filtered10.txt", quote = F, row.names = F)

plotA <- ggplot(working_ac_heart, aes(x = sample_id, y = log10_ref)) + geom_violin() + ggtitle("Heart") + ylab("log10 Ref Count")
plotB <- ggplot(working_ac_heart, aes(x = sample_id, y = log10_alt)) + geom_violin() + ggtitle("Heart") + ylab("log10 Alt Count")

plotA + plotB


######################################
#### Filtering the allele counts #####
######################################

#working_ac_heart_filtered <- working_ac_heart %>%
#  filter(
#    (sample_id == "F1") |
#      (sample_id == "F2") |
#      (sample_id == "F3") |
#      (sample_id == "M1"))
#If I do this, I lose a lot of data but it has the most accurate representation of the data 
#since the other males have much more extreme values for the allele counts for reference: alternate allele counts


working_ac_heart %>%
  subset(sample_id == "F1") %>% pull(ref_count) %>% quantile(probs = 0.2) # 15
working_ac_heart %>%
  subset(sample_id == "F1") %>% pull(ref_count) %>% quantile(probs = 0.8) # 83
working_ac_heart %>%
  subset(sample_id == "F1") %>% pull(alt_count) %>% quantile(probs = 0.2) # 12
working_ac_heart %>%
  subset(sample_id == "F1") %>% pull(alt_count) %>% quantile(probs = 0.8) # 42

working_ac_heart %>%
  subset(sample_id == "F2") %>% pull(ref_count) %>% quantile(probs = 0.2) # 15
working_ac_heart %>%
  subset(sample_id == "F2") %>% pull(ref_count) %>% quantile(probs = 0.8) # 77
working_ac_heart %>%
  subset(sample_id == "F2") %>% pull(alt_count) %>% quantile(probs = 0.2) # 12
working_ac_heart %>%
  subset(sample_id == "F2") %>% pull(alt_count) %>% quantile(probs = 0.8) # 39

working_ac_heart %>%
  subset(sample_id == "F3") %>% pull(ref_count) %>% quantile(probs = 0.2) # 15
working_ac_heart %>%
  subset(sample_id == "F3") %>% pull(ref_count) %>% quantile(probs = 0.8) # 84
working_ac_heart %>%
  subset(sample_id == "F3") %>% pull(alt_count) %>% quantile(probs = 0.2) # 13
working_ac_heart %>%
  subset(sample_id == "F3") %>% pull(alt_count) %>% quantile(probs = 0.8) # 44


working_ac_heart %>%
  subset(sample_id == "M1") %>% pull(ref_count) %>% quantile(probs = 0.2) # 15
working_ac_heart %>%
  subset(sample_id == "M1") %>% pull(ref_count) %>% quantile(probs = 0.8) # 85
working_ac_heart %>%
  subset(sample_id == "M1") %>% pull(alt_count) %>% quantile(probs = 0.2) # 12
working_ac_heart %>%
  subset(sample_id == "M1") %>% pull(alt_count) %>% quantile(probs = 0.8) # 43

working_ac_heart %>%
  subset(sample_id == "M2") %>% pull(ref_count) %>% quantile(probs = 0.2) # 19
working_ac_heart %>%
  subset(sample_id == "M2") %>% pull(ref_count) %>% quantile(probs = 0.8) # 167
working_ac_heart %>%
  subset(sample_id == "M2") %>% pull(alt_count) %>% quantile(probs = 0.2) # 13
working_ac_heart %>%
  subset(sample_id == "M2") %>% pull(alt_count) %>% quantile(probs = 0.8) # 68

working_ac_heart %>%
  subset(sample_id == "M3") %>% pull(ref_count) %>% quantile(probs = 0.2) # 18
working_ac_heart %>%
  subset(sample_id == "M3") %>% pull(ref_count) %>% quantile(probs = 0.8) # 144
working_ac_heart %>%
  subset(sample_id == "M3") %>% pull(alt_count) %>% quantile(probs = 0.2) # 13
working_ac_heart %>%
  subset(sample_id == "M3") %>% pull(alt_count) %>% quantile(probs = 0.8) # 64

###


working_ac_heart_filtered <- working_ac_heart %>%
  filter(
    (sample_id == "F1" & ref_count < 42 & alt_count > 12) |
    (sample_id == "F2" & ref_count < 42 & alt_count > 12) |
    (sample_id == "F3" & ref_count < 42 & alt_count > 12) |
    (sample_id == "M1" & ref_count < 41 & alt_count > 12) |
    (sample_id == "M2" & ref_count < 41 & alt_count > 12) |
    (sample_id == "M3" & ref_count < 41 & alt_count > 12))

plotC <- ggplot(working_ac_heart_filtered, aes(x = sample_id, y = log10_ref)) + geom_violin() + ggtitle("Heart Filtered") + ylab("log10 Ref Count")
plotD <- ggplot(working_ac_heart_filtered, aes(x = sample_id, y = log10_alt)) + geom_violin() + ggtitle("Heart Filtered") + ylab("log10 Alt Count")

plotC + plotD

#working_ac_heart_filtered <- working_ac_heart

#####################################

#Attach chrom_type:
gene_list <- c("PDE8B", "unc5da", "ENSPREG00000019378", "ENSPREG00000019435", "ENSPREG00000019478", "ENSPREG00000019540",
               "ENSPREG00000019560", "ENSPREG00000019673", "st8sia3", "ENSPREG00000019838", "WDR7", "ENSPREG00000019881",
               "ccng2", "ccni", "ENSPREG00000019970", "ENSPREG00000019988", "ENSPREG00000020020", "ENSPREG00000020032",
               "ENSPREG00000020049", "ENSPREG00000020126", "ENSPREG00000020131", "TMC2", "zfand5b", "ccnb1", "CACNA1B",
               "vldlr", "ENSPREG00000001200", "c7b", "ENSPREG00000001446", "ENSPREG00000001498", "ENSPREG00000001901",
               "ENSPREG00000005507", "ENSPREG00000005501", "ENSPREG00000005616", "ENSPREG00000005759", "stoml2", "PIGO",
               "ENSPREG00000006339", "ENSPREG00000006387", "ENSPREG00000006501", 
               "KCNV2", "ENSPREG00000019622", "BRCC3", "dnajc25", "NARS1")

### Everything past "ENSPREG00000006501" were genes found from the previous HI_SNPs pipeline with the INDELS removed
### And it was all filtered through snp_ilter_by_individ_CLEANED_FINAL.R that weren't already included
# ENSPREG00000019446 = KCNV2
# ENSPREG00000019635 = BRCC3
# ENSPREG00000019738 = dnajc25
# ENSPREG00000020098 = NARS1

add_chrom_column <- function(dfA, dfB) {
  matched_idx <- match(dfA$geneid, dfB$GENE_NAME)
  dfA$chrom <- dfB$CHROM[matched_idx]
  return(dfA)
}

working_ac_heart_filtered <- add_chrom_column(working_ac_heart_filtered, gene_info)
working_ac_heart_filtered$chrom_type <- ifelse(working_ac_heart_filtered$chrom == "LG12", "Sex Chromosome", "Autosome")
working_ac_heart_filtered$chrom_type[working_ac_heart_filtered$geneid %in% gene_list] <- "Impacted Gene"

working_ac_heart_A <- subset(working_ac_heart_filtered, chrom_type == "Autosome")
working_ac_heart_sc <- subset(working_ac_heart_filtered, chrom_type == "Sex Chromosome")
working_ac_heart_HI <- subset(working_ac_heart_filtered, chrom_type == "Impacted Gene")

combo_working_ac_heart <- rbind(working_ac_heart_A, working_ac_heart_HI, working_ac_heart_sc)

celltype_ase <- combo_working_ac_heart %>%
  group_by(sample_id, geneid) %>%
  summarise(ref_cluster = mean(ref_count), alt_cluster = mean(alt_count))
celltype_ase$sex <- ifelse(grepl("^F", celltype_ase$sample_id), "F", "M")
celltype_ase$total_counts <- celltype_ase$ref_cluster + celltype_ase$alt_cluster
celltype_ase$mar <- pmax(celltype_ase$ref_cluster, celltype_ase$alt_cluster)/(celltype_ase$total_counts)

#celltype_ase_unfilt <- add_chrom_column(celltype_ase, gene_info)
#celltype_ase_unfilt$chrom_type <- ifelse(celltype_ase_unfilt$chrom == "LG12", "Sex Chromosome", "Autosome")
#celltype_ase_unfilt$chrom_type[celltype_ase_unfilt$geneid %in% gene_list] <- "Impacted Gene"

#celltype_ase_IG <- subset(celltype_ase_unfilt, chrom_type == "Impacted Gene")
#celltype_ase_filt <- celltype_ase_unfilt %>%
#  filter(chrom_type != "Impacted Gene") %>%
#  group_by(geneid) %>% filter(n_distinct(sample_id[sex == "F"]) >= 3, n_distinct(sample_id[sex == "M"]) >= 3,) %>%
#  ungroup()

#celltype_ase <- rbind(celltype_ase_IG, celltype_ase_filt)

plot_ase <- celltype_ase %>%
  group_by(geneid, sex) %>%
  summarise(avg_mar = mean(mar))
plot_ase <- plot_ase %>%
  arrange(desc(avg_mar), desc(geneid))

plot_ase <- add_chrom_column(plot_ase, gene_info)
plot_ase$chrom_type <- ifelse(plot_ase$chrom == "LG12", "Sex Chromosome", "Autosome")
plot_ase$chrom_type[plot_ase$geneid %in% gene_list] <- "Impacted Gene"

heart_impactgenes <- subset(plot_ase, chrom_type == "Impacted Gene")
heart_autosome <- subset(plot_ase, chrom_type == "Autosome")
heart_sexchromo <- subset(plot_ase, chrom_type == "Sex Chromosome")

### Plotting:

heart_impactgenes %>% filter(sex == "F") %>% pull(avg_mar) %>% median()
heart_impactgenes %>% filter(sex == "M") %>% pull(avg_mar) %>% median()
heart_autosome %>% filter(sex == "F") %>% pull(avg_mar) %>% median()
heart_autosome %>% filter(sex == "M") %>% pull(avg_mar) %>% median()
heart_sexchromo %>% filter(sex == "F") %>% pull(avg_mar) %>% median()
heart_sexchromo %>% filter(sex == "M") %>% pull(avg_mar) %>% median()

plot_heart_IG <- ggplot(heart_impactgenes, aes(x = avg_mar)) + geom_density(aes(fill = sex), alpha = 0.5) + 
  scale_fill_manual(values = c("firebrick3", "cornflowerblue")) +
  #scale_fill_manual(values = c("firebrick3", "hotpink", "pink", "navy", "blue", "cornflowerblue")) + 
  labs(fill = "Sex") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 12)) +
  ylab("Density") + xlab("Major Allele Ratio") + ggtitle("Heart (Impacted Genes)") +
  geom_vline(aes(xintercept=0.56),color="firebrick4",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5539019),color="navy",linetype="dashed",size=0.8) +
  theme_classic() + theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    plot.title = element_text(size = 20, hjust = 0.5)
  )
plot_heart_IG

plot_heart_A <- ggplot(heart_autosome, aes(x = avg_mar)) + geom_density(aes(fill = sex), alpha = 0.5) + 
  scale_fill_manual(values = c("firebrick3", "cornflowerblue")) +
  #scale_fill_manual(values = c("firebrick3", "hotpink", "pink", "navy", "blue", "cornflowerblue")) + 
  labs(fill = "Sex") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 12)) +
  ylab("Density") + xlab("Major Allele Ratio") + ggtitle("Heart (Autosomes)") +
  geom_vline(aes(xintercept=0.5610803),color="firebrick4",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.569438),color="navy",linetype="dashed",size=0.8) +
  theme_classic() + theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    plot.title = element_text(size = 20, hjust = 0.5)
  )
plot_heart_A

plot_heart_SC <- ggplot(heart_sexchromo, aes(x = avg_mar)) + geom_density(aes(fill = sex), alpha = 0.5) + 
  scale_fill_manual(values = c("firebrick3", "cornflowerblue")) +
  #scale_fill_manual(values = c("firebrick3", "hotpink", "pink", "navy", "blue", "cornflowerblue")) + 
  labs(fill = "Sex") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 12)) +
  ylab("Density") + xlab("Major Allele Ratio") + ggtitle("Heart (Sex Chromosome)") +
  geom_vline(aes(xintercept=0.5735248),color="firebrick4",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5675631),color="navy",linetype="dashed",size=0.8) +
  theme_classic() + theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    plot.title = element_text(size = 20, hjust = 0.5)
  )
plot_heart_SC

heart_plots <- (plot_heart_IG | plot_heart_A | plot_heart_SC) + plot_layout(guides = "collect")
heart_plots

wilcox.test(heart_impactgenes$avg_mar ~ heart_impactgenes$sex, p.adjust.methods = "bonferroni")
wilcox.test(heart_autosome$avg_mar ~ heart_autosome$sex, p.adjust.methods = "bonferroni")
wilcox.test(heart_sexchromo$avg_mar ~ heart_sexchromo$sex, p.adjust.methods = "bonferroni")

#############################################################################
########### Run DDS on all celltypes before splitting them up: ##############
#############################################################################


# Convert to long format
ref_df_temp <- as.data.frame(ref_matrix_all) %>%
  tibble::rownames_to_column("gene_id") %>%
  pivot_longer(-gene_id, names_to = "cell", values_to = "ref")

alt_df_temp <- as.data.frame(alt_matrix_all) %>%
  tibble::rownames_to_column("gene_id") %>%
  pivot_longer(-gene_id, names_to = "cell", values_to = "alt")

colnames(ref_df_temp)[1] <- "cell"
colnames(ref_df_temp)[2] <- "gene_id"
colnames(alt_df_temp)[1] <- "cell"
colnames(alt_df_temp)[2] <- "gene_id"

# Merge counts + assign sex
combined_df_temp <- ref_df_temp %>%
  mutate(alt = alt_df_temp$alt,
         sex = ifelse(grepl("F", cell), "F", "M")) 

# Reshape to long format
long_df_temp <- combined_df_temp %>%
  pivot_longer(cols = c(ref, alt), names_to = "allele", values_to = "count") %>%
  mutate(
    cell = factor(cell),
    allele = factor(allele, levels = c("ref", "alt")),
    sex = factor(sex)
  ) %>%
  unite(sample, cell, allele, sep = "_", remove = FALSE) %>%
  na.omit()

male_df <- filter(long_df_temp, sex == "M")
female_df <- filter(long_df_temp, sex == "F")

# Run dds - male:
count_mat_male <- male_df %>%
  select(gene_id, sample, count) %>%
  pivot_wider(names_from = sample, values_from = count, values_fill = 0) %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

coldata_male <- male_df %>%
  distinct(sample, allele, cell) %>%
  as.data.frame()
rownames(coldata_male) <- coldata_male$sample

dds_male <- DESeqDataSetFromMatrix(
  countData = count_mat_male,
  colData = coldata_male,
  design = ~ allele
)

dds_male <-estimateSizeFactors(dds_male, type = "poscounts")
dds_male <- DESeq(dds_male)

res_male <- results(dds_male, contrast = c("allele", "alt", "ref"))
res_df_male <- as.data.frame(res_male) %>%
  filter(log2FoldChange != 0)
sig_res_df_male <- as.data.frame(res_male) %>%
  filter(log2FoldChange != 0, padj < 0.05)
write.csv(res_df_male, "scASE_results/heart_male_all.csv", row.names = T)
write.csv(sig_res_df_male, "scASE_results/heart_male_all_sig.csv", row.names = T)

# Run dds - female:
count_mat_female <- female_df %>%
  select(gene_id, sample, count) %>%
  pivot_wider(names_from = sample, values_from = count, values_fill = 0) %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

coldata_female <- female_df %>%
  distinct(sample, allele, cell) %>%
  as.data.frame()
rownames(coldata_female) <- coldata_female$sample

dds_female <- DESeqDataSetFromMatrix(
  countData = count_mat_female,
  colData = coldata_female,
  design = ~ allele
)

dds_female <-estimateSizeFactors(dds_female, type = "poscounts")
dds_female <- DESeq(dds_female)

res_female <- results(dds_female, contrast = c("allele", "alt", "ref"))
res_df_female <- as.data.frame(res_female) %>%
  filter(log2FoldChange != 0)
sig_res_df_female <- as.data.frame(res_female) %>%
  filter(log2FoldChange != 0, padj < 0.05)
write.csv(res_df_female, "scASE_results/heart_female_all.csv", row.names = T)
write.csv(sig_res_df_female, "scASE_results/heart_female_all_sig.csv", row.names = T)


#######################################################################
########### To run this for all celltypes in your matrix ##############
#######################################################################

# Write a function so you can pull out every matrix from your ^_pb list:
# Separate out male vs. female

run_scASE_by_sex <- function(cell_type, ref_list, alt_list, output_dir = "scASE_results/heart/") {
  message("Processing: ", cell_type)
  
  # Extract matrices
  ref_mat <- as.matrix(ref_list[[cell_type]])
  alt_mat <- as.matrix(alt_list[[cell_type]])
  
  # Convert to long format
  ref_df <- as.data.frame(ref_mat) %>%
    tibble::rownames_to_column("gene_id") %>%
    pivot_longer(-gene_id, names_to = "cell", values_to = "ref")
  
  alt_df <- as.data.frame(alt_mat) %>%
    tibble::rownames_to_column("gene_id") %>%
    pivot_longer(-gene_id, names_to = "cell", values_to = "alt")
  
  # Merge counts + assign sex
  combined_df <- ref_df %>%
    mutate(alt = alt_df$alt,
           sex = ifelse(grepl("F", cell), "F", "M")) 
  
  # Reshape to long format
  long_df <- combined_df %>%
    pivot_longer(cols = c(ref, alt), names_to = "allele", values_to = "count") %>%
    mutate(
      cell = factor(cell),
      allele = factor(allele, levels = c("ref", "alt")),
      sex = factor(sex)
    ) %>%
    unite(sample, cell, allele, sep = "_", remove = FALSE) %>%
    na.omit()
  
  # Define helper for each sex subset
  run_sex_specific_DESeq2 <- function(sex_label) {
    sex_df <- filter(long_df, sex == sex_label)
    
    count_mat <- sex_df %>%
      select(gene_id, sample, count) %>%
      pivot_wider(names_from = sample, values_from = count, values_fill = 0) %>%
      column_to_rownames("gene_id") %>%
      as.matrix()
    
    coldata <- sex_df %>%
      distinct(sample, allele, cell) %>%
      as.data.frame()
    rownames(coldata) <- coldata$sample
    
    dds <- DESeqDataSetFromMatrix(
      countData = count_mat,
      colData = coldata,
      design = ~ allele
    )
    
    dds <- DESeq(dds)
    res <- results(dds, contrast = c("allele", "alt", "ref"))
    res_df <- as.data.frame(res) %>%
      filter(log2FoldChange != 0)
    sig_res_df <- as.data.frame(res) %>%
      filter(log2FoldChange != 0, padj < 0.05)
    
    # Save
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    write.csv(res_df, file.path(output_dir, paste0(cell_type, "_", sex_label, ".csv")), row.names = T)
    write.csv(sig_res_df, file.path(output_dir, paste0(cell_type, "_sig_", sex_label, ".csv")), row.names = T)
    
    return(res_df)
  }
  
  list(
    F = run_sex_specific_DESeq2("F"),
    M = run_sex_specific_DESeq2("M")
  )
}

cell_types <- names(ref_pb_heart)

results_by_sex <- lapply(cell_types, function(ct) {
  run_scASE_by_sex(ct, ref_pb_heart, alt_pb_heart)
})
names(results_by_sex) <- cell_types


#########################################
############ To visualize: ##############
#########################################

### Gene list that has HI SNPs

gene_list <- c("PDE8B", "unc5da", "ENSPREG00000019378", "ENSPREG00000019435", "ENSPREG00000019478", "ENSPREG00000019540",
               "ENSPREG00000019560", "ENSPREG00000019673", "st8sia3", "ENSPREG00000019838", "WDR7", "ENSPREG00000019881",
               "ccng2", "ccni", "ENSPREG00000019970", "ENSPREG00000019988", "ENSPREG00000020020", "ENSPREG00000020032",
               "ENSPREG00000020049", "ENSPREG00000020126", "ENSPREG00000020131", "TMC2", "zfand5b", "ccnb1", "CACNA1B",
               "vldlr", "ENSPREG00000001200", "c7b", "ENSPREG00000001446", "ENSPREG00000001498", "ENSPREG00000001901",
               "ENSPREG00000005507", "ENSPREG00000005501", "ENSPREG00000005616", "ENSPREG00000005759", "stoml2", "PIGO",
               "ENSPREG00000006339", "ENSPREG00000006387", "ENSPREG00000006501", 
               "KCNV2", "ENSPREG00000019622", "BRCC3", "dnajc25", "NARS1")

### Everything past "ENSPREG00000006501" were genes found from the previous HI_SNPs pipeline with the INDELS removed
### And it was all filtered through snp_ilter_by_individ_CLEANED_FINAL.R that weren't already included
# ENSPREG00000019446 = KCNV2
# ENSPREG00000019635 = BRCC3
# ENSPREG00000019738 = dnajc25
# ENSPREG00000020098 = NARS1

add_chrom_column <- function(dfA, dfB) {
  matched_idx <- match(dfA$X, dfB$GENE_NAME)
  dfA$chrom <- dfB$CHROM[matched_idx]
  return(dfA)
}


    # 1. Cardiomyocyte

card_F <- read.csv("scASE_results/heart/Cardiomyocyte_F.csv")
card_M <- read.csv("scASE_results/heart/Cardiomyocyte_M.csv")
card_F$sex = "F"
card_M$sex = "M"
card_combo <- rbind(card_F, card_M)

card_combo$prop <- (2^card_combo$log2FoldChange)/(1+2^card_combo$log2FoldChange)
card_combo$mar <- pmax(card_combo$prop, 1 - card_combo$prop)
card_combo <- add_chrom_column(card_combo, gene_info)
card_combo$chrom_type <- ifelse(card_combo$chrom == "LG12", "Sex Chromosome", "Autosome")
card_combo$chrom_type[card_combo$X %in% gene_list] <- "Impacted Gene"
card_combo_M <- subset(card_combo, sex == "M")
card_combo_F <- subset(card_combo, sex == "F")

card_combo_M %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.6082877
card_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.6216497
card_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.5967469


card_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.5812739
card_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.5281052
card_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.5970551

write.csv(card_combo_M, "scASE_results/heart/split_chromo/cardiomyocyte_M.csv", quote = FALSE, row.names = FALSE)
write.csv(card_combo_F, "scASE_results/heart/split_chromo/cardiomyocyte_F.csv", quote = FALSE, row.names = FALSE)


plot_card_M <- ggplot(card_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 8)) +
  geom_vline(aes(xintercept=0.6082877),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6216497),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5967469),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Heart (Cardiomyocyte M)")
plot_card_M

plot_card_F <- ggplot(card_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 3) + 
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 8)) +
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  geom_vline(aes(xintercept=0.5812739),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5281052),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5970551),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Heart (Cardiomyocyte F)")

plot_card_F


    # 2. Endothelial

endo_F <- read.csv("scASE_results/heart/Endothelial_F.csv")
endo_M <- read.csv("scASE_results/heart/Endothelial_M.csv")
endo_F$sex = "F"
endo_M$sex = "M"
endo_combo <- rbind(endo_F, endo_M)

endo_combo$prop <- (2^endo_combo$log2FoldChange)/(1+2^endo_combo$log2FoldChange)
endo_combo$mar <- pmax(endo_combo$prop, 1 - endo_combo$prop)
endo_combo <- add_chrom_column(endo_combo, gene_info)
endo_combo$chrom_type <- ifelse(endo_combo$chrom == "LG12", "Sex Chromosome", "Autosome")
endo_combo$chrom_type[endo_combo$X %in% gene_list] <- "Impacted Gene"
endo_combo_M <- subset(endo_combo, sex == "M")
endo_combo_F <- subset(endo_combo, sex == "F")

endo_combo_M %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.6220719
endo_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.5746337
endo_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6047359


endo_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.5795079
endo_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.5994241
endo_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.5911904

write.csv(endo_combo_M, "scASE_results/heart/split_chromo/Endothelial_M.csv", quote = FALSE, row.names = FALSE)
write.csv(endo_combo_F, "scASE_results/heart/split_chromo/Endothelial_F.csv", quote = FALSE, row.names = FALSE)


plot_endo_M <- ggplot(endo_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 4) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 8)) +
  geom_vline(aes(xintercept=0.6220719),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5746337),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6047359),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Heart (Endothelial M)")
plot_endo_M

plot_endo_F <- ggplot(endo_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 2) + 
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 8)) +
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  geom_vline(aes(xintercept=0.5795079),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5994241),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5911904),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Heart (Endothelial F)")

plot_endo_F


    # 3. Erythrocyte

eryth_F <- read.csv("scASE_results/heart/Erythrocyte_F.csv")
eryth_M <- read.csv("scASE_results/heart/Erythrocyte_M.csv")
eryth_F$sex = "F"
eryth_M$sex = "M"
eryth_combo <- rbind(eryth_F, eryth_M)

eryth_combo$prop <- (2^eryth_combo$log2FoldChange)/(1+2^eryth_combo$log2FoldChange)
eryth_combo$mar <- pmax(eryth_combo$prop, 1 - eryth_combo$prop)
eryth_combo <- add_chrom_column(eryth_combo, gene_info)
eryth_combo$chrom_type <- ifelse(eryth_combo$chrom == "LG12", "Sex Chromosome", "Autosome")
eryth_combo$chrom_type[eryth_combo$X %in% gene_list] <- "Impacted Gene"
eryth_combo_M <- subset(eryth_combo, sex == "M")
eryth_combo_F <- subset(eryth_combo, sex == "F")

eryth_combo_M %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.6024734
eryth_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.6181506
eryth_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6103918


eryth_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.5833821
eryth_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.5462837
eryth_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.5920632

write.csv(eryth_combo_M, "scASE_results/heart/split_chromo/Erythrocyte_M.csv", quote = FALSE, row.names = FALSE)
write.csv(eryth_combo_F, "scASE_results/heart/split_chromo/Erythrocyte_F.csv", quote = FALSE, row.names = FALSE)


plot_eryth_M <- ggplot(eryth_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 8)) +
  geom_vline(aes(xintercept=0.6024734),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6181506),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6103918),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Heart (Erythrocyte M)")
plot_eryth_M

plot_eryth_F <- ggplot(eryth_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.5) + 
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 8)) +
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  geom_vline(aes(xintercept=0.5833821),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5462837),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5920632),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Heart (Erythrocyte F)")

plot_eryth_F


    # 4. Fibroblast

fibro_F <- read.csv("scASE_results/heart/Fibroblast_F.csv")
fibro_M <- read.csv("scASE_results/heart/Fibroblast_M.csv")
fibro_F$sex = "F"
fibro_M$sex = "M"
fibro_combo <- rbind(fibro_F, fibro_M)

fibro_combo$prop <- (2^fibro_combo$log2FoldChange)/(1+2^fibro_combo$log2FoldChange)
fibro_combo$mar <- pmax(fibro_combo$prop, 1 - fibro_combo$prop)
fibro_combo <- add_chrom_column(fibro_combo, gene_info)
fibro_combo$chrom_type <- ifelse(fibro_combo$chrom == "LG12", "Sex Chromosome", "Autosome")
fibro_combo$chrom_type[fibro_combo$X %in% gene_list] <- "Impacted Gene"
fibro_combo_M <- subset(fibro_combo, sex == "M")
fibro_combo_F <- subset(fibro_combo, sex == "F")

fibro_combo_M %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.5995243
fibro_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.5832634
fibro_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6025659


fibro_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.5693826
fibro_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.5683947
fibro_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.5864793

write.csv(fibro_combo_M, "scASE_results/heart/split_chromo/Fibroblast_M.csv", quote = FALSE, row.names = FALSE)
write.csv(fibro_combo_F, "scASE_results/heart/split_chromo/Fibroblast_F.csv", quote = FALSE, row.names = FALSE)


plot_fibro_M <- ggplot(fibro_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 8)) +
  geom_vline(aes(xintercept=0.5995243),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5832634),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6025659),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Heart (Fibroblast M)")
plot_fibro_M

plot_fibro_F <- ggplot(fibro_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.5) + 
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 8)) +
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  geom_vline(aes(xintercept=0.5693826),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5683947),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5864793),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Heart (Fibroblast F)")

plot_fibro_F


    # 5. Granulocyte

granu_F <- read.csv("scASE_results/heart/Granulocyte_F.csv")
granu_M <- read.csv("scASE_results/heart/Granulocyte_M.csv")
granu_F$sex = "F"
granu_M$sex = "M"
granu_combo <- rbind(granu_F, granu_M)

granu_combo$prop <- (2^granu_combo$log2FoldChange)/(1+2^granu_combo$log2FoldChange)
granu_combo$mar <- pmax(granu_combo$prop, 1 - granu_combo$prop)
granu_combo <- add_chrom_column(granu_combo, gene_info)
granu_combo$chrom_type <- ifelse(granu_combo$chrom == "LG12", "Sex Chromosome", "Autosome")
granu_combo$chrom_type[granu_combo$X %in% gene_list] <- "Impacted Gene"
granu_combo_M <- subset(granu_combo, sex == "M")
granu_combo_F <- subset(granu_combo, sex == "F")

granu_combo_M %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.6256981
granu_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.5720113
granu_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6504745

granu_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) #  0.6373487
granu_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.6101953
granu_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.620629

write.csv(granu_combo_M, "scASE_results/heart/split_chromo/Granulocyte_M.csv", quote = FALSE, row.names = FALSE)
write.csv(granu_combo_F, "scASE_results/heart/split_chromo/Granulocyte_F.csv", quote = FALSE, row.names = FALSE)


plot_granu_M <- ggplot(granu_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 8)) +
  geom_vline(aes(xintercept=0.6256981),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5720113),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6504745),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Heart (Granulocyte M)")
plot_granu_M

plot_granu_F <- ggplot(granu_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 2) + 
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 8)) +
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  geom_vline(aes(xintercept=0.6373487),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6101953),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.620629),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Heart (Granulocyte F)")
plot_granu_F


    # 6. Macrophage

macro_F <- read.csv("scASE_results/heart/Macrophage_F.csv")
macro_M <- read.csv("scASE_results/heart/Macrophage_M.csv")
macro_F$sex = "F"
macro_M$sex = "M"
macro_combo <- rbind(macro_F, macro_M)

macro_combo$prop <- (2^macro_combo$log2FoldChange)/(1+2^macro_combo$log2FoldChange)
macro_combo$mar <- pmax(macro_combo$prop, 1 - macro_combo$prop)
macro_combo <- add_chrom_column(macro_combo, gene_info)
macro_combo$chrom_type <- ifelse(macro_combo$chrom == "LG12", "Sex Chromosome", "Autosome")
macro_combo$chrom_type[macro_combo$X %in% gene_list] <- "Impacted Gene"
macro_combo_M <- subset(macro_combo, sex == "M")
macro_combo_F <- subset(macro_combo, sex == "F")

macro_combo_M %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.6024155
macro_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.6008689
macro_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6253291

macro_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) #  0.5729485
macro_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.5596384
macro_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.5776933

write.csv(macro_combo_M, "scASE_results/heart/split_chromo/Macrophage_M.csv", quote = FALSE, row.names = FALSE)
write.csv(macro_combo_F, "scASE_results/heart/split_chromo/Macrophage_F.csv", quote = FALSE, row.names = FALSE)


plot_macro_M <- ggplot(macro_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 8)) +
  geom_vline(aes(xintercept=0.6024155),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6008689),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6253291),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Heart (Macrophage M)")
plot_macro_M

plot_macro_F <- ggplot(macro_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.5) + 
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 8)) +
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  geom_vline(aes(xintercept=0.5729485),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5596384),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5776933),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Heart (Macrophage F)")
plot_macro_F


    # 7. Immune (Other Immune)

immune_F <- read.csv("scASE_results/heart/Other immune_F.csv")
immune_M <- read.csv("scASE_results/heart/Other immune_M.csv")
immune_F$sex = "F"
immune_M$sex = "M"
immune_combo <- rbind(immune_F, immune_M)

immune_combo$prop <- (2^immune_combo$log2FoldChange)/(1+2^immune_combo$log2FoldChange)
immune_combo$mar <- pmax(immune_combo$prop, 1 - immune_combo$prop)
immune_combo <- add_chrom_column(immune_combo, gene_info)
immune_combo$chrom_type <- ifelse(immune_combo$chrom == "LG12", "Sex Chromosome", "Autosome")
immune_combo$chrom_type[immune_combo$X %in% gene_list] <- "Impacted Gene"
immune_combo_M <- subset(immune_combo, sex == "M")
immune_combo_F <- subset(immune_combo, sex == "F")

immune_combo_M %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.6337553
immune_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.5368422
immune_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) #  0.6334726

immune_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) #  0.5928954
immune_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.558279
immune_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6087169

write.csv(immune_combo_M, "scASE_results/heart/split_chromo/other_immune_M.csv", quote = FALSE, row.names = FALSE)
write.csv(immune_combo_F, "scASE_results/heart/split_chromo/other_immune_F.csv", quote = FALSE, row.names = FALSE)


plot_immune_M <- ggplot(immune_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 8)) +
  geom_vline(aes(xintercept=0.6337553),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5368422),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept= 0.6334726),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Heart (Other immune cells M)")
plot_immune_M

plot_immune_F <- ggplot(immune_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.5) + 
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 8)) +
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  geom_vline(aes(xintercept=0.5928954),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.558279),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6087169),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Heart (Other immune cells F)")
plot_immune_F


    # 8. Smooth Muscle Cell

smooth_F <- read.csv("scASE_results/heart/Smooth muscle cell_F.csv")
smooth_M <- read.csv("scASE_results/heart/Smooth muscle cell_M.csv")
smooth_F$sex = "F"
smooth_M$sex = "M"
smooth_combo <- rbind(smooth_F, smooth_M)

smooth_combo$prop <- (2^smooth_combo$log2FoldChange)/(1+2^smooth_combo$log2FoldChange)
smooth_combo$mar <- pmax(smooth_combo$prop, 1 - smooth_combo$prop)
smooth_combo <- add_chrom_column(smooth_combo, gene_info)
smooth_combo$chrom_type <- ifelse(smooth_combo$chrom == "LG12", "Sex Chromosome", "Autosome")
smooth_combo$chrom_type[smooth_combo$X %in% gene_list] <- "Impacted Gene"
smooth_combo_M <- subset(smooth_combo, sex == "M")
smooth_combo_F <- subset(smooth_combo, sex == "F")

smooth_combo_M %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.6219939
smooth_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.6647288
smooth_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) #  0.6311738

smooth_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) #  0.6092227
smooth_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.6136686
smooth_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6084712

write.csv(smooth_combo_M, "scASE_results/heart/split_chromo/Smooth_muscle_M.csv", quote = FALSE, row.names = FALSE)
write.csv(smooth_combo_F, "scASE_results/heart/split_chromo/Smooth_muscle_F.csv", quote = FALSE, row.names = FALSE)


plot_smooth_M <- ggplot(smooth_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.2) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 8)) +
  geom_vline(aes(xintercept=0.6219939),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6647288),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept= 0.6311738),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Heart (Smooth muscle cells M)")
plot_smooth_M

plot_smooth_F <- ggplot(smooth_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.3) + 
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 8)) +
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  geom_vline(aes(xintercept=0.6092227),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6136686),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6084712),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Heart (Smooth muscle cells F)")
plot_smooth_F


    # 9. T lymophocyte

tlymph_F <- read.csv("scASE_results/heart/T lymphocyte_F.csv")
tlymph_M <- read.csv("scASE_results/heart/T lymphocyte_M.csv")
tlymph_F$sex = "F"
tlymph_M$sex = "M"
tlymph_combo <- rbind(tlymph_F, tlymph_M)

tlymph_combo$prop <- (2^tlymph_combo$log2FoldChange)/(1+2^tlymph_combo$log2FoldChange)
tlymph_combo$mar <- pmax(tlymph_combo$prop, 1 - tlymph_combo$prop)
tlymph_combo <- add_chrom_column(tlymph_combo, gene_info)
tlymph_combo$chrom_type <- ifelse(tlymph_combo$chrom == "LG12", "Sex Chromosome", "Autosome")
tlymph_combo$chrom_type[tlymph_combo$X %in% gene_list] <- "Impacted Gene"
tlymph_combo_M <- subset(tlymph_combo, sex == "M")
tlymph_combo_F <- subset(tlymph_combo, sex == "F")

tlymph_combo_M %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.615145
tlymph_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.5725352
tlymph_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) #  0.6271312

tlymph_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) #  0.5704611
tlymph_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) #  0.5534633
tlymph_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.5746974

write.csv(tlymph_combo_M, "scASE_results/heart/split_chromo/T_lymphocyte_M.csv", quote = FALSE, row.names = FALSE)
write.csv(tlymph_combo_F, "scASE_results/heart/split_chromo/T_lymphocyte_F.csv", quote = FALSE, row.names = FALSE)


plot_tlymph_M <- ggplot(tlymph_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.2) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 8)) +
  geom_vline(aes(xintercept=0.615145),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5725352),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept= 0.6271312),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Heart (T lymphocyte M)")
plot_tlymph_M

plot_tlymph_F <- ggplot(tlymph_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 2.5) + 
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 8)) +
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  geom_vline(aes(xintercept=0.5704611),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5534633),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5746974),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Heart (T lymphocyte F)")
plot_tlymph_F



#### Combo your plots ####

library(patchwork)

heart_celltype_plots_M <- wrap_plots(list(plot_card_M, plot_endo_M, plot_eryth_M, plot_fibro_M, plot_granu_M, 
                                         plot_immune_M, plot_macro_M, plot_smooth_M, plot_tlymph_M), ncol = 3, nrow = 3) +
  plot_layout(guides = "collect")

heart_celltype_plots_M

heart_celltype_plots_F <- wrap_plots(list(plot_card_F, plot_endo_F, plot_eryth_F, plot_fibro_F, plot_granu_F, 
                                          plot_immune_F, plot_macro_F, plot_smooth_F, plot_tlymph_F), ncol = 3, nrow = 3) +
  plot_layout(guides = "collect")

heart_celltype_plots_F



######################################################
######## Combine celltypes to test for ASE ###########
######################################################


all_celltype_M <- rbind(card_combo_M, endo_combo_M, eryth_combo_M, fibro_combo_M, granu_combo_M, immune_combo_M,
                        macro_combo_M, smooth_combo_M, tlymph_combo_M)

all_celltype_F <- rbind(card_combo_F, endo_combo_F, eryth_combo_F, fibro_combo_F, granu_combo_F, immune_combo_F,
                        macro_combo_F, smooth_combo_F, tlymph_combo_F)


all_celltype_M_means <- all_celltype_M %>%
  group_by(X, chrom_type) %>%
  summarise(mar_mean = mean(mar), .groups = "drop")

all_celltype_F_means <- all_celltype_F %>%
  group_by(X, chrom_type) %>%
  summarise(mar_mean = mean(mar), .groups = "drop")



all_celltype_M_means %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar_mean)) #  0.623
all_celltype_M_means %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar_mean)) # 0.616
all_celltype_M_means %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar_mean)) # 0.622

all_celltype_F_means %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar_mean)) # 0.590
all_celltype_F_means %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar_mean)) # 0.568
all_celltype_F_means %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar_mean)) # 0.590


plot_celltype_M <- ggplot(all_celltype_M_means, aes(x = mar_mean)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.2) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 15)) +
  geom_vline(aes(xintercept=0.623),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.616),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept= 0.622),color="darkorchid4",linetype="dashed",size=0.8) +
  theme_classic() + ggtitle("Heart M (all celltypes)")
plot_celltype_M


plot_celltype_F <- ggplot(all_celltype_F_means, aes(x = mar_mean)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.5) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 15)) +
  geom_vline(aes(xintercept=0.590),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.568),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept= 0.590),color="darkorchid4",linetype="dashed",size=0.8) +
  theme_classic() + ggtitle("Heart F (all celltypes)")
plot_celltype_F


combo_everything <- wrap_plots(plot_celltype_M, plot_celltype_F) +
  plot_layout(guides = "collect")

combo_everything



