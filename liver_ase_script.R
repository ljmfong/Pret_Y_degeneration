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

#### 1. Load in required packages ####

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
##### 2. Load in Iulia's matrices ####

load("iulia_scRNA/liver_umap_3_rename_clusters_dbrmv.RData")

#### 3. Get differential gene expression for separate cell types ####
counts_liver <- liver_umap_3_rename_clusters_dbrmv@assays$RNA@counts
metadata_liver <- liver_umap_3_rename_clusters_dbrmv@meta.data
metadata_liver$CellType <- factor(liver_umap_3_rename_clusters_dbrmv@active.ident)

########

#allelecount_Lib2 <- read.table("scrna_scAlleleCount/complete_matrice_for_R/Lib2_counts.txt", sep = " ", header = T)
#allelecount_Lib12 <- read.table("scrna_scAlleleCount/complete_matrice_for_R/Lib12_counts.txt", sep = " ", header = T)
#allelecount_Lib14 <- read.table("scrna_scAlleleCount/complete_matrice_for_R/Lib14_counts.txt", sep = " ", header = T)
#allelecount_Lib8 <- read.table("scrna_scAlleleCount/complete_matrice_for_R/Lib8_counts.txt", sep = " ", header = T)
#allelecount_Lib9 <- read.table("scrna_scAlleleCount/complete_matrice_for_R/Lib9_counts.txt", sep = " ", header = T)
#allelecount_Lib10 <- read.table("scrna_scAlleleCount/complete_matrice_for_R/Lib10_counts.txt", sep = " ", header = T)

allelecount_Lib2 <- read.table("cellSNP_scAlleleCount/complete_matrice_for_R/Lib2_counts.txt", sep = " ", header = T)
allelecount_Lib12 <- read.table("cellSNP_scAlleleCount/complete_matrice_for_R/Lib12_counts.txt", sep = " ", header = T)
allelecount_Lib14 <- read.table("cellSNP_scAlleleCount/complete_matrice_for_R/Lib14_counts.txt", sep = " ", header = T)
allelecount_Lib8 <- read.table("cellSNP_scAlleleCount/complete_matrice_for_R/Lib8_counts.txt", sep = " ", header = T)
allelecount_Lib9 <- read.table("cellSNP_scAlleleCount/complete_matrice_for_R/Lib9_counts.txt", sep = " ", header = T)
allelecount_Lib10 <- read.table("cellSNP_scAlleleCount/complete_matrice_for_R/Lib10_counts.txt", sep = " ", header = T)

merged_liver_AC <- rbind(allelecount_Lib2, allelecount_Lib12, allelecount_Lib14, 
                         allelecount_Lib8, allelecount_Lib9, allelecount_Lib10)
colnames(merged_liver_AC)[3] <- "cells"


#### Attach gene info and find matching gene_name for positional information:

gene_info <- read.table("new_Ret_SNPs/gene_boundary_w_gene_names.bed", header = T, sep = "\t")


attached_gene_AC <- merged_liver_AC %>%
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
#write.table(cleaned_df, "scASE_results/liver_cleaned_df_gene_attached.txt",  quote = F, row.names = F)
write.table(cleaned_df, "cellSNP_scAlleleCount/R_results/liver_cleaned_df_gene_attached.txt",  quote = F, row.names = F)
#if loading in:
#cleaned_df <- read.table(file = "scASE_results/liver_cleaned_df_gene_attached.txt", header = T)
cleaned_df <- read.table(file = "cellSNP_scAlleleCount/R_results/liver_cleaned_df_gene_attached.txt", header = T)

######### Attach to your SingleCellExperiment

#OG sce_liver command:
sce_liver <- SingleCellExperiment(assays=list(counts=counts_liver), colData=metadata_liver)

# Determine the number of clusters (nk) and the cluster names (kids), and the number of samples (ns) and sample names (sids)
colData(sce_liver)$sample_id <- as.factor(colData(sce_liver)$sample)
colData(sce_liver)$sex_id <- as.factor(colData(sce_liver)$sex)
colData(sce_liver)$cluster_id <- as.factor(colData(sce_liver)$CellType)
kids_liver <- purrr::set_names(levels(sce_liver$cluster_id))
nk_liver <- length(kids_liver)
sids_liver <- purrr::set_names(levels(sce_liver$sample_id))
ns_liver <- length(sids_liver)

# Turn named vector into a numeric vector
n_cells_liver <- as.numeric(table(sce_liver$sample_id))

# Reorder samples (rows) of the metadata to match the order of the sample names
m_liver <- match(sids_liver, sce_liver$sample_id)

# Create sample level metadata by combining the reordered metadata with the number of cells
ei_liver <- data.frame(colData(sce_liver)[m_liver, ], n_cells_liver, row.names=NULL) %>% select(-"cluster_id")

# Remove lowly expressed genes
sce_liver <- sce_liver[rowSums(counts(sce_liver)>1)>=10,]	

# Aggregate counts per sample_id and cluster_id
groups_liver <- colData(sce_liver)[, c("cluster_id", "sample_id")]
pb_liver <- aggregate.Matrix(t(counts(sce_liver)), groupings=groups_liver, fun="sum")


##############


summed_ac_counts <- cleaned_df %>%
  group_by(cells, geneid) %>%
  summarise(ref_count = sum(ref, na.rm = T), alt_count = sum(alt, na.rm = T), .groups = "drop")

any(duplicated(summed_ac_counts$cells)) #This should be TRUE


summed_ac_counts <- summed_ac_counts %>%
  mutate(cells = str_replace(cells, "^Lib8", "F1_Lib8"),
         cells = str_replace(cells, "^Lib9", "F2_Lib9"),
         cells = str_replace(cells, "^Lib10", "F3_Lib10"),
         cells = str_replace(cells, "^Lib2", "M1_Lib2"),
         cells = str_replace(cells, "^Lib12", "M2_Lib12"),
         cells = str_replace(cells, "^Lib14", "M3_Lib14"))

### Cluster your summed_ac_counts into cell-types
#All the info for the match is in groups_livers

groups_liver_name <- data.frame(groups_liver) 
groups_liver_name <- groups_liver_name %>% rownames_to_column(var = "cells")

match(groups_liver_name$cells, summed_ac_counts$cells)

group_name_ac_liver <- summed_ac_counts %>%
  left_join(groups_liver_name, by = "cells")

cleaned_group_name_ac_liver <- na.omit(group_name_ac_liver)


################ Turn into matrix ####################
############## save as dgCMatrix #####################


# Build ONE big sparse matrix (cells x genes) for ref_count & alt_count
ref_matrix_all <- cleaned_group_name_ac_liver %>%
  select(cells, geneid, ref_count) %>%
  group_by(cells, geneid) %>%
  summarise(count = sum(ref_count), .groups = "drop") %>%
  pivot_wider(names_from = geneid, values_from = count, values_fill = 0) %>%
  column_to_rownames("cells") %>%
  as.matrix() %>%
  Matrix(sparse = TRUE)


alt_matrix_all <- cleaned_group_name_ac_liver %>%
  select(cells, geneid, alt_count) %>%
  group_by(cells, geneid) %>%
  summarise(count = sum(alt_count), .groups = "drop") %>%
  pivot_wider(names_from = geneid, values_from = count, values_fill = 0) %>%
  column_to_rownames("cells") %>%
  as.matrix() %>%
  Matrix(sparse = TRUE)

cell_metadata <- cleaned_group_name_ac_liver %>%
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

ref_splitf_liver <- sapply(stringr::str_split(rownames(aggregated_ref), pattern="_", n=2), '[', 1) # Check that each cluster is present in all samples
alt_splitf_liver <- sapply(stringr::str_split(rownames(aggregated_alt), pattern="_", n=2), '[', 1) # Check that each cluster is present in all samples

ref_pb_liver <- split.data.frame(aggregated_ref, factor(ref_splitf_liver)) %>% lapply(function(u) set_colnames(t(u), stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))
alt_pb_liver <- split.data.frame(aggregated_alt, factor(alt_splitf_liver)) %>% lapply(function(u) set_colnames(t(u), stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))


##############################################################################
########### Find MAR on all celltypes before splitting them up: ##############
##############################################################################

working_ac_liver <- summed_ac_counts
working_ac_liver <- working_ac_liver %>%
  filter(ref_count >= 10 & alt_count >= 10)
working_ac_liver$log10_ref <- log10(working_ac_liver$ref_count)
working_ac_liver$log10_alt <- log10(working_ac_liver$alt_count)

working_ac_liver <- working_ac_liver %>%
  mutate(sample_id = case_when( 
    grepl("^M1", cells) ~ "M1",
    grepl("^M2", cells) ~ "M2",
    grepl("^M3", cells) ~ "M3",
    grepl("^F1", cells) ~ "F1",
    grepl("^F2", cells) ~ "F2",
    grepl("^F3", cells) ~ "F3"))

#write.table(working_ac_liver, "scASE_results/liver_filtered10.txt", quote = F, row.names = F)
write.table(working_ac_liver, "cellSNP_scAlleleCount/R_results/liver_filtered10.txt", quote = F, row.names = F)

plotA <- ggplot(working_ac_liver, aes(x = sample_id, y = log10_ref)) + geom_violin() + ggtitle("Liver") + ylab("log10 Ref Count")
plotB <- ggplot(working_ac_liver, aes(x = sample_id, y = log10_alt)) + geom_violin() + ggtitle("Liver") + ylab("log10 Alt Count")

plotA + plotB

#Distribution of the allele counts look equal across males vs. female samples

######################################
#### Filtering the allele counts #####
######################################

working_ac_liver %>%
  subset(sample_id == "F1") %>% pull(ref_count) %>% quantile(probs = 0.2) # 15
working_ac_liver %>%
  subset(sample_id == "F1") %>% pull(ref_count) %>% quantile(probs = 0.8) # 85
working_ac_liver %>%
  subset(sample_id == "F1") %>% pull(alt_count) %>% quantile(probs = 0.2) # 12
working_ac_liver %>%
  subset(sample_id == "F1") %>% pull(alt_count) %>% quantile(probs = 0.8) # 43

working_ac_liver %>%
  subset(sample_id == "F2") %>% pull(ref_count) %>% quantile(probs = 0.2) # 16
working_ac_liver %>%
  subset(sample_id == "F2") %>% pull(ref_count) %>% quantile(probs = 0.8) # 93
working_ac_liver %>%
  subset(sample_id == "F2") %>% pull(alt_count) %>% quantile(probs = 0.2) # 12
working_ac_liver %>%
  subset(sample_id == "F2") %>% pull(alt_count) %>% quantile(probs = 0.8) # 47

working_ac_liver %>%
  subset(sample_id == "F3") %>% pull(ref_count) %>% quantile(probs = 0.2) # 15
working_ac_liver %>%
  subset(sample_id == "F3") %>% pull(ref_count) %>% quantile(probs = 0.8) # 87
working_ac_liver %>%
  subset(sample_id == "F3") %>% pull(alt_count) %>% quantile(probs = 0.2) # 12
working_ac_liver %>%
  subset(sample_id == "F3") %>% pull(alt_count) %>% quantile(probs = 0.8) # 43


working_ac_liver %>%
  subset(sample_id == "M1") %>% pull(ref_count) %>% quantile(probs = 0.2) # 15
working_ac_liver %>%
  subset(sample_id == "M1") %>% pull(ref_count) %>% quantile(probs = 0.8) # 81
working_ac_liver %>%
  subset(sample_id == "M1") %>% pull(alt_count) %>% quantile(probs = 0.2) # 12
working_ac_liver %>%
  subset(sample_id == "M1") %>% pull(alt_count) %>% quantile(probs = 0.8) # 41

working_ac_liver %>%
  subset(sample_id == "M2") %>% pull(ref_count) %>% quantile(probs = 0.2) # 16
working_ac_liver %>%
  subset(sample_id == "M2") %>% pull(ref_count) %>% quantile(probs = 0.8) # 106
working_ac_liver %>%
  subset(sample_id == "M2") %>% pull(alt_count) %>% quantile(probs = 0.2) # 13
working_ac_liver %>%
  subset(sample_id == "M2") %>% pull(alt_count) %>% quantile(probs = 0.8) # 48
 
working_ac_liver %>%
  subset(sample_id == "M3") %>% pull(ref_count) %>% quantile(probs = 0.2) # 16
working_ac_liver %>%
  subset(sample_id == "M3") %>% pull(ref_count) %>% quantile(probs = 0.8) # 105
working_ac_liver %>%
  subset(sample_id == "M3") %>% pull(alt_count) %>% quantile(probs = 0.2) # 13
working_ac_liver %>%
  subset(sample_id == "M3") %>% pull(alt_count) %>% quantile(probs = 0.8) # 51


#Pick the lowest distribution for males and females separately:

working_ac_liver_filtered <- working_ac_liver %>%
  filter(
    (sample_id == "F1" & ref_count < 43 & alt_count > 12) |
    (sample_id == "F2" & ref_count < 43 & alt_count > 12) |
    (sample_id == "F3" & ref_count < 43 & alt_count > 12) |
    (sample_id == "M1" & ref_count < 41 & alt_count > 12) |
    (sample_id == "M2" & ref_count < 41 & alt_count > 12) |
    (sample_id == "M3" & ref_count < 41 & alt_count > 12))

plotC <- ggplot(working_ac_liver_filtered, aes(x = sample_id, y = log10_ref)) + geom_violin(trim = F) + ggtitle("liver Filtered") + ylab("log10 Ref Count")
plotD <- ggplot(working_ac_liver_filtered, aes(x = sample_id, y = log10_alt)) + geom_violin(trim = F) + ggtitle("liver Filtered") + ylab("log10 Alt Count")

plotC + plotD

#working_ac_liver_filtered <- working_ac_liver

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

working_ac_liver_filtered <- add_chrom_column(working_ac_liver_filtered, gene_info)
working_ac_liver_filtered$chrom_type <- ifelse(working_ac_liver_filtered$chrom == "LG12", "Sex Chromosome", "Autosome")
working_ac_liver_filtered$chrom_type[working_ac_liver_filtered$geneid %in% gene_list] <- "Impacted Gene"

working_ac_liver_A <- subset(working_ac_liver_filtered, chrom_type == "Autosome")
working_ac_liver_sc <- subset(working_ac_liver_filtered, chrom_type == "Sex Chromosome")
working_ac_liver_HI <- subset(working_ac_liver_filtered, chrom_type == "Impacted Gene")

combo_working_ac_liver <- rbind(working_ac_liver_A, working_ac_liver_HI, working_ac_liver_sc)


celltype_ase <- combo_working_ac_liver %>%
  group_by(sample_id, geneid) %>%
  summarise(ref_cluster = mean(ref_count), alt_cluster = mean(alt_count))
celltype_ase$sex <- ifelse(grepl("^F", celltype_ase$sample_id), "F", "M")
celltype_ase$total_counts <- celltype_ase$ref_cluster + celltype_ase$alt_cluster
celltype_ase$mar <- pmax(celltype_ase$ref_cluster, celltype_ase$alt_cluster)/(celltype_ase$total_counts)

plot_ase <- celltype_ase %>%
  group_by(geneid, sex) %>%
  summarise(avg_mar = mean(mar))
plot_ase <- plot_ase %>%
  arrange(desc(avg_mar), desc(geneid))

plot_ase <- add_chrom_column(plot_ase, gene_info)
plot_ase$chrom_type <- ifelse(plot_ase$chrom == "LG12", "Sex Chromosome", "Autosome")
plot_ase$chrom_type[plot_ase$geneid %in% gene_list] <- "Impacted Gene"

liver_impactgenes <- subset(plot_ase, chrom_type == "Impacted Gene")
liver_autosome <- subset(plot_ase, chrom_type == "Autosome")
liver_sexchromo <- subset(plot_ase, chrom_type == "Sex Chromosome")

### Plotting:

liver_impactgenes %>% filter(sex == "F") %>% pull(avg_mar) %>% median()
liver_impactgenes %>% filter(sex == "M") %>% pull(avg_mar) %>% median()
liver_autosome %>% filter(sex == "F") %>% pull(avg_mar) %>% median()
liver_autosome %>% filter(sex == "M") %>% pull(avg_mar) %>% median()
liver_sexchromo %>% filter(sex == "F") %>% pull(avg_mar) %>% median()
liver_sexchromo %>% filter(sex == "M") %>% pull(avg_mar) %>% median()

plot_liver_IG <- ggplot(liver_impactgenes, aes(x = avg_mar)) + geom_density(aes(fill = sex), alpha = 0.5) + 
  scale_fill_manual(values = c("firebrick3", "cornflowerblue")) +
  #scale_fill_manual(values = c("firebrick3", "hotpink", "pink", "navy", "blue", "cornflowerblue")) + 
  labs(fill = "Sex") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 12)) +
  ylab("Density") + xlab("Major Allele Ratio") + ggtitle("Liver (Impacted Genes)") +
  geom_vline(aes(xintercept=0.5774847),color="firebrick4",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5445291),color="navy",linetype="dashed",size=0.8) +
  theme_classic() + theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    plot.title = element_text(size = 20, hjust = 0.5)
  )
plot_liver_IG

plot_liver_A <- ggplot(liver_autosome, aes(x = avg_mar)) + geom_density(aes(fill = sex), alpha = 0.5) + 
  scale_fill_manual(values = c("firebrick3", "cornflowerblue")) +
  #scale_fill_manual(values = c("firebrick3", "hotpink", "pink", "navy", "blue", "cornflowerblue")) + 
  labs(fill = "Sex") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 12)) +
  ylab("Density") + xlab("Major Allele Ratio") + ggtitle("Liver (Autosomes)") +
  geom_vline(aes(xintercept=0.558442),color="firebrick4",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5589841),color="navy",linetype="dashed",size=0.8) +
  theme_classic() + theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    plot.title = element_text(size = 20, hjust = 0.5)
  )
plot_liver_A

plot_liver_SC <- ggplot(liver_sexchromo, aes(x = avg_mar)) + geom_density(aes(fill = sex), alpha = 0.5, adjust = 1.1) + 
  scale_fill_manual(values = c("firebrick3", "cornflowerblue")) +
  #scale_fill_manual(values = c("firebrick3", "hotpink", "pink", "navy", "blue", "cornflowerblue")) + 
  labs(fill = "Sex") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 12)) +
  ylab("Density") + xlab("Major Allele Ratio") + ggtitle("Liver (Sex Chromosome)") +
  geom_vline(aes(xintercept=0.5673451),color="firebrick4",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.56),color="navy",linetype="dashed",size=0.8) +
  theme_classic() + theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    plot.title = element_text(size = 20, hjust = 0.5)
  )
plot_liver_SC

liver_plots <- (plot_liver_IG | plot_liver_A | plot_liver_SC) + plot_layout(guides = "collect")
liver_plots

wilcox.test(liver_impactgenes$avg_mar ~ liver_impactgenes$sex, p.adjust.methods = "bonferroni")
wilcox.test(liver_autosome$avg_mar ~ liver_autosome$sex, p.adjust.methods = "bonferroni")
wilcox.test(liver_sexchromo$avg_mar ~ liver_sexchromo$sex, p.adjust.methods = "bonferroni")


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
write.csv(res_df_male, "scASE_results/liver_male_all.csv", row.names = T)
write.csv(sig_res_df_male, "scASE_results/liver_male_all_sig.csv", row.names = T)

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
write.csv(res_df_female, "scASE_results/liver_female_all.csv", row.names = T)
write.csv(sig_res_df_female, "scASE_results/liver_female_all_sig.csv", row.names = T)

########### To run this for all celltypes in your matrix ##############

# Write a function so you can pull out every matrix from your ^_pb list:
# Separate out male vs. female

run_scASE_by_sex <- function(cell_type, ref_list, alt_list, output_dir = "scASE_results/liver/") {
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

cell_types <- names(ref_pb_liver)

results_by_sex <- lapply(cell_types, function(ct) {
  run_scASE_by_sex(ct, ref_pb_liver, alt_pb_liver)
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

    # 1. B lymphocyte

blymph_F <- read.csv("scASE_results/liver/B lymphocyte_F.csv")
blymph_M <- read.csv("scASE_results/liver/B lymphocyte_M.csv")
blymph_F$sex = "F"
blymph_M$sex = "M"
blymph_combo <- rbind(blymph_F, blymph_M)

blymph_combo$prop <- (2^blymph_combo$log2FoldChange)/(1+2^blymph_combo$log2FoldChange)
blymph_combo$mar <- pmax(blymph_combo$prop, 1 - blymph_combo$prop)
blymph_combo <- add_chrom_column(blymph_combo, gene_info)
blymph_combo$chrom_type <- ifelse(blymph_combo$chrom == "LG12", "Sex Chromosome", "Autosome")
blymph_combo$chrom_type[blymph_combo$X %in% gene_list] <- "Impacted Gene"
blymph_combo_M <- subset(blymph_combo, sex == "M")
blymph_combo_F <- subset(blymph_combo, sex == "F")

blymph_combo_M %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) #  0.6019171
blymph_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.5657691
blymph_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6099494


blymph_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.5924706
blymph_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.5836793
blymph_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6104933

write.csv(blymph_combo_M, "scASE_results/liver/split_chromo/blymphiomyocyte_M.csv", quote = FALSE, row.names = FALSE)
write.csv(blymph_combo_F, "scASE_results/liver/split_chromo/blymphiomyocyte_F.csv", quote = FALSE, row.names = FALSE)


plot_blymph_M <- ggplot(blymph_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.5) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 10)) +
  geom_vline(aes(xintercept= 0.6019171),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5657691),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6099494),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Liver (B lymphocyte M)")
plot_blymph_M

plot_blymph_F <- ggplot(blymph_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 2.2) + 
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 10)) +
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  geom_vline(aes(xintercept=0.5924706),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5836793),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6104933),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Liver (B lymphocyte F)")
plot_blymph_F


    # 2. Biliary epithelial

Bepi_F <- read.csv("scASE_results/liver/Biliary epithelial_F.csv")
Bepi_M <- read.csv("scASE_results/liver/Biliary epithelial_M.csv")
Bepi_F$sex = "F"
Bepi_M$sex = "M"
Bepi_combo <- rbind(Bepi_F, Bepi_M)

Bepi_combo$prop <- (2^Bepi_combo$log2FoldChange)/(1+2^Bepi_combo$log2FoldChange)
Bepi_combo$mar <- pmax(Bepi_combo$prop, 1 - Bepi_combo$prop)
Bepi_combo <- add_chrom_column(Bepi_combo, gene_info)
Bepi_combo$chrom_type <- ifelse(Bepi_combo$chrom == "LG12", "Sex Chromosome", "Autosome")
Bepi_combo$chrom_type[Bepi_combo$X %in% gene_list] <- "Impacted Gene"
Bepi_combo_M <- subset(Bepi_combo, sex == "M")
Bepi_combo_F <- subset(Bepi_combo, sex == "F")

Bepi_combo_M %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) #  0.6073038
Bepi_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.6189324
Bepi_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6248117


Bepi_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.6035871
Bepi_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.5843997
Bepi_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6033362

write.csv(Bepi_combo_M, "scASE_results/liver/split_chromo/Biliary_epithelial_M.csv", quote = FALSE, row.names = FALSE)
write.csv(Bepi_combo_F, "scASE_results/liver/split_chromo/Biliary_epithelial_F.csv", quote = FALSE, row.names = FALSE)


plot_Bepi_M <- ggplot(Bepi_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.5) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 10)) +
  geom_vline(aes(xintercept= 0.6073038),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6189324),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6248117),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Liver (Biliary epithelial M)")
plot_Bepi_M

plot_Bepi_F <- ggplot(Bepi_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.5) + 
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 10)) +
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  geom_vline(aes(xintercept=0.6035871),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5843997),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6033362),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Liver (Biliary epithelial F)")
plot_Bepi_F


    # 3. Endothelial

endo_F <- read.csv("scASE_results/liver/Endothelial_F.csv")
endo_M <- read.csv("scASE_results/liver/Endothelial_M.csv")
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
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) #  0.657261
endo_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.5353918
endo_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6629474


endo_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.6085748
endo_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.6184344
endo_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.625322

write.csv(endo_combo_M, "scASE_results/liver/split_chromo/Endothelial_M.csv", quote = FALSE, row.names = FALSE)
write.csv(endo_combo_F, "scASE_results/liver/split_chromo/Endothelial_F.csv", quote = FALSE, row.names = FALSE)


plot_endo_M <- ggplot(endo_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.2) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 10)) +
  geom_vline(aes(xintercept= 0.657261),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5353918),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6629474),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Liver (Endothelial M)")
plot_endo_M

plot_endo_F <- ggplot(endo_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.2) + 
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 10)) +
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  geom_vline(aes(xintercept=0.6085748),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6184344),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.625322),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Liver (Endothelial F)")
plot_endo_F

    # 4. Erythrocyte

eryth_F <- read.csv("scASE_results/liver/Erythrocyte_F.csv")
eryth_M <- read.csv("scASE_results/liver/Erythrocyte_M.csv")
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
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) #  0.586531
eryth_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.5847883
eryth_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.5997831


eryth_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.5804182
eryth_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.5168376
eryth_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.5734478

write.csv(eryth_combo_M, "scASE_results/liver/split_chromo/Erythrocyte_M.csv", quote = FALSE, row.names = FALSE)
write.csv(eryth_combo_F, "scASE_results/liver/split_chromo/Erythrocyte_F.csv", quote = FALSE, row.names = FALSE)


plot_eryth_M <- ggplot(eryth_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 4) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 10)) +
  geom_vline(aes(xintercept= 0.586531),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5847883),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5997831),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Liver (Erythrocyte M)")
plot_eryth_M

plot_eryth_F <- ggplot(eryth_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 10) + 
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 10)) +
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  geom_vline(aes(xintercept=0.5804182),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5168376),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5734478),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Liver (Erythrocyte F)")
plot_eryth_F

    # 5. Hepatocyte 1

hep1_F <- read.csv("scASE_results/liver/Hepatocyte 1_F.csv")
hep1_M <- read.csv("scASE_results/liver/Hepatocyte 1_M.csv")
hep1_F$sex = "F"
hep1_M$sex = "M"
hep1_combo <- rbind(hep1_F, hep1_M)

hep1_combo$prop <- (2^hep1_combo$log2FoldChange)/(1+2^hep1_combo$log2FoldChange)
hep1_combo$mar <- pmax(hep1_combo$prop, 1 - hep1_combo$prop)
hep1_combo <- add_chrom_column(hep1_combo, gene_info)
hep1_combo$chrom_type <- ifelse(hep1_combo$chrom == "LG12", "Sex Chromosome", "Autosome")
hep1_combo$chrom_type[hep1_combo$X %in% gene_list] <- "Impacted Gene"
hep1_combo_M <- subset(hep1_combo, sex == "M")
hep1_combo_F <- subset(hep1_combo, sex == "F")

hep1_combo_M %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) #  0.5602128
hep1_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.5581981
hep1_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.5687892


hep1_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.5585791
hep1_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.541677
hep1_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.5648684

write.csv(hep1_combo_M, "scASE_results/liver/split_chromo/Hep_1_M.csv", quote = FALSE, row.names = FALSE)
write.csv(hep1_combo_F, "scASE_results/liver/split_chromo/Hep_1_F.csv", quote = FALSE, row.names = FALSE)


plot_hep1_M <- ggplot(hep1_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.5) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 10)) +
  geom_vline(aes(xintercept= 0.5602128),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5687892),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5997831),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Liver (Hepatocyte 1 M)")
plot_hep1_M

plot_hep1_F <- ggplot(hep1_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 3) + 
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 10)) +
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  geom_vline(aes(xintercept=0.5585791),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.541677),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5648684),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Liver (Hepatocyte 1 F)")
plot_hep1_F


    # 6. Hepatocyte 2

hep2_F <- read.csv("scASE_results/liver/Hepatocyte 2_F.csv")
hep2_M <- read.csv("scASE_results/liver/Hepatocyte 2_M.csv")
hep2_F$sex = "F"
hep2_M$sex = "M"
hep2_combo <- rbind(hep2_F, hep2_M)

hep2_combo$prop <- (2^hep2_combo$log2FoldChange)/(1+2^hep2_combo$log2FoldChange)
hep2_combo$mar <- pmax(hep2_combo$prop, 1 - hep2_combo$prop)
hep2_combo <- add_chrom_column(hep2_combo, gene_info)
hep2_combo$chrom_type <- ifelse(hep2_combo$chrom == "LG12", "Sex Chromosome", "Autosome")
hep2_combo$chrom_type[hep2_combo$X %in% gene_list] <- "Impacted Gene"
hep2_combo_M <- subset(hep2_combo, sex == "M")
hep2_combo_F <- subset(hep2_combo, sex == "F")

hep2_combo_M %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) #  0.5926708
hep2_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.5979701
hep2_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.5941548


hep2_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.6058033
hep2_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.6322297
hep2_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6280799

write.csv(hep2_combo_M, "scASE_results/liver/split_chromo/Hep_2_M.csv", quote = FALSE, row.names = FALSE)
write.csv(hep2_combo_F, "scASE_results/liver/split_chromo/Hep_2_F.csv", quote = FALSE, row.names = FALSE)


plot_hep2_M <- ggplot(hep2_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 10)) +
  geom_vline(aes(xintercept= 0.5926708),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5979701),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5997831),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Liver (Hepatocyte 2 M)")
plot_hep2_M

plot_hep2_F <- ggplot(hep2_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.2) + 
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 10)) +
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  geom_vline(aes(xintercept=0.6058033),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6322297),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6280799),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Liver (Hepatocyte 2 F)")
plot_hep2_F


    # 7. Macrophage

macro_F <- read.csv("scASE_results/liver/Macrophage_F.csv")
macro_M <- read.csv("scASE_results/liver/Macrophage_M.csv")
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
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) #  0.5861138
macro_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.5477191
macro_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.604809


macro_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.58322
macro_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.5944026
macro_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.5843036

write.csv(macro_combo_M, "scASE_results/liver/split_chromo/Macrophage_M.csv", quote = FALSE, row.names = FALSE)
write.csv(macro_combo_F, "scASE_results/liver/split_chromo/Macrophage_F.csv", quote = FALSE, row.names = FALSE)


plot_macro_M <- ggplot(macro_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 10)) +
  geom_vline(aes(xintercept= 0.5861138),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5477191),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.604809),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Liver (Macrophage M)")
plot_macro_M

plot_macro_F <- ggplot(macro_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.2) + 
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 10)) +
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  geom_vline(aes(xintercept=0.58322),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5944026),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5843036),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Liver (Macrophage F)")
plot_macro_F

    # 8. Other Immune

immune_F <- read.csv("scASE_results/liver/Other immune_F.csv")
immune_M <- read.csv("scASE_results/liver/Other immune_M.csv")
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
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) #  0.5848302
immune_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.6120273
immune_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.5812354


immune_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.5680566
immune_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.5729093
immune_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.5710537

write.csv(immune_combo_M, "scASE_results/liver/split_chromo/immune_M.csv", quote = FALSE, row.names = FALSE)
write.csv(immune_combo_F, "scASE_results/liver/split_chromo/immune_F.csv", quote = FALSE, row.names = FALSE)


plot_immune_M <- ggplot(immune_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 10)) +
  geom_vline(aes(xintercept= 0.5848302),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6120273),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5812354),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Liver (Other immune cells M)")
plot_immune_M

plot_immune_F <- ggplot(immune_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.2) + 
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 10)) +
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  geom_vline(aes(xintercept=0.5680566),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5729093),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5710537),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Liver (Other immune cells F)")
plot_immune_F

    # 9. T lymphocyte

tlymph_F <- read.csv("scASE_results/liver/T lymphocyte_F.csv")
tlymph_M <- read.csv("scASE_results/liver/T lymphocyte_M.csv")
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
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) #  0.5971311
tlymph_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.5948722
tlymph_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6049997


tlymph_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.6020717
tlymph_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.6047785
tlymph_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6133788

write.csv(tlymph_combo_M, "scASE_results/liver/split_chromo/Tlymph_M.csv", quote = FALSE, row.names = FALSE)
write.csv(tlymph_combo_F, "scASE_results/liver/split_chromo/Tlymph_F.csv", quote = FALSE, row.names = FALSE)


plot_tlymph_M <- ggplot(tlymph_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 10)) +
  geom_vline(aes(xintercept= 0.5971311),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5948722),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6049997),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Liver (T lymphocyte M)")
plot_tlymph_M

plot_tlymph_F <- ggplot(tlymph_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.5) + 
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 10)) +
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  geom_vline(aes(xintercept=0.6020717),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6047785),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6133788),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Liver (T lymphocyte F)")
plot_tlymph_F

#############################################

#### Combo your plots ####

library(patchwork)

liver_celltype_plots_M <- wrap_plots(list(plot_blymph_M, plot_Bepi_M, plot_endo_M, plot_eryth_M, plot_hep1_M, 
                                          plot_hep2_M, plot_immune_M, plot_macro_M, plot_tlymph_M), ncol = 3, nrow = 3) +
  plot_layout(guides = "collect")

liver_celltype_plots_M

liver_celltype_plots_F <- wrap_plots(list(plot_blymph_F, plot_Bepi_F, plot_endo_F, plot_eryth_F, plot_hep1_F, 
                                          plot_hep2_F, plot_immune_F, plot_macro_F, plot_tlymph_F), ncol = 3, nrow = 3) +
  plot_layout(guides = "collect")

liver_celltype_plots_F



######################################################
######## Combine celltypes to test for ASE ###########
######################################################


all_celltype_M <- rbind(blymph_combo_M, Bepi_combo_M, eryth_combo_M, endo_combo_M, hep1_combo_M, immune_combo_M,
                        macro_combo_M, hep2_combo_M, tlymph_combo_M)

all_celltype_F <- rbind(blymph_combo_F, Bepi_combo_F, eryth_combo_F, endo_combo_F, hep1_combo_F, immune_combo_F,
                        macro_combo_F, hep2_combo_F, tlymph_combo_F)


all_celltype_M_means <- all_celltype_M %>%
  group_by(X, chrom_type) %>%
  summarise(mar_mean = mean(mar), .groups = "drop")

all_celltype_F_means <- all_celltype_F %>%
  group_by(X, chrom_type) %>%
  summarise(mar_mean = mean(mar), .groups = "drop")



all_celltype_M_means %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar_mean)) #  0.604
all_celltype_M_means %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar_mean)) # 0.598
all_celltype_M_means %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar_mean)) # 0.609

all_celltype_F_means %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar_mean)) # 0.595
all_celltype_F_means %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar_mean)) # 0.605
all_celltype_F_means %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar_mean)) # 0.595


plot_celltype_M <- ggplot(all_celltype_M_means, aes(x = mar_mean)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.5) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 15)) +
  geom_vline(aes(xintercept=0.604),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.598),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept= 0.609),color="darkorchid4",linetype="dashed",size=0.8) +
  theme_classic() + ggtitle("Liver M (all celltypes)")
plot_celltype_M


plot_celltype_F <- ggplot(all_celltype_F_means, aes(x = mar_mean)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.5) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 15)) +
  geom_vline(aes(xintercept=0.595),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.605),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept= 0.595),color="darkorchid4",linetype="dashed",size=0.8) +
  theme_classic() + ggtitle("Liver F (all celltypes)")
plot_celltype_F


combo_everything <- wrap_plots(plot_celltype_M, plot_celltype_F) +
  plot_layout(guides = "collect")

combo_everything

