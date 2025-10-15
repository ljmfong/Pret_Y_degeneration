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

load("iulia_scRNA/skin_umap_3_rename_clusters_dbrmv.RData")

counts_skin <- skin_umap_3_rename_clusters_dbrmv@assays$RNA@counts
metadata_skin <- skin_umap_3_rename_clusters_dbrmv@meta.data
metadata_skin$CellType <- factor(skin_umap_3_rename_clusters_dbrmv@active.ident)

#allelecount_Lib19 <- read.table("scrna_scAlleleCount/complete_matrice_for_R/Lib19_counts.txt", sep = " ", header = T)
#allelecount_Lib18 <- read.table("scrna_scAlleleCount/complete_matrice_for_R/Lib18_counts.txt", sep = " ", header = T)
#allelecount_Lib26 <- read.table("scrna_scAlleleCount/complete_matrice_for_R/Lib26_counts.txt", sep = " ", header = T)
#allelecount_Lib32 <- read.table("scrna_scAlleleCount/complete_matrice_for_R/Lib32_counts.txt", sep = " ", header = T)
#allelecount_Lib33 <- read.table("scrna_scAlleleCount/complete_matrice_for_R/Lib33_counts.txt", sep = " ", header = T)
#allelecount_Lib34 <- read.table("scrna_scAlleleCount/complete_matrice_for_R/Lib34_counts.txt", sep = " ", header = T)

allelecount_Lib19 <- read.table("cellSNP_scAlleleCount/complete_matrice_for_R/Lib19_counts.txt", sep = " ", header = T)
allelecount_Lib18 <- read.table("cellSNP_scAlleleCount/complete_matrice_for_R/Lib18_counts.txt", sep = " ", header = T)
allelecount_Lib26 <- read.table("cellSNP_scAlleleCount/complete_matrice_for_R/Lib26_counts.txt", sep = " ", header = T)
allelecount_Lib32 <- read.table("cellSNP_scAlleleCount/complete_matrice_for_R/Lib32_counts.txt", sep = " ", header = T)
allelecount_Lib33 <- read.table("cellSNP_scAlleleCount/complete_matrice_for_R/Lib33_counts.txt", sep = " ", header = T)
allelecount_Lib34 <- read.table("cellSNP_scAlleleCount/complete_matrice_for_R/Lib34_counts.txt", sep = " ", header = T)

merged_skin_AC <- rbind(allelecount_Lib19, allelecount_Lib18, allelecount_Lib26, 
                        allelecount_Lib32, allelecount_Lib33, allelecount_Lib34)
colnames(merged_skin_AC)[3] <- "cells"


#### Attach gene info and find matching gene_name for positional information:

gene_info <- read.table("new_Ret_SNPs/gene_boundary_w_gene_names.bed", header = T, sep = "\t")

attached_gene_AC <- merged_skin_AC %>%
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
#write.table(cleaned_df, "scASE_results/skin_cleaned_df_gene_attached.txt", quote = F, row.names = F)
write.table(cleaned_df, "cellSNP_scAlleleCount/R_results/skin_cleaned_df_gene_attached.txt", quote = F, row.names = F)
#If loading in:
#cleaned_df <- read.table(file = "scASE_results/skin_cleaned_df_gene_attached.txt", header = T)
cleaned_df <- read.table(file = "cellSNP_scAlleleCount/R_results/skin_cleaned_df_gene_attached.txt", header = T)

#matched_indices <- match(metadata_skin$cells, attached_gene_AC$cell_name)
#matched_df <- attached_gene_AC[matched_indices, ]


#full_df <- metadata_skin %>%
#  select(cells) %>% 
#  left_join(cleaned_df, by = c("cells" = "cell_name")) %>%
#  mutate(
#    ref = replace_na(ref, 0),
#    alt = replace_na(alt, 0),
#    CHROM = replace_na(CHROM, "0"),
#    POS = replace_na(POS, 0),
#    geneid = replace_na(geneid, "0")
#  )


######### Attach to your SingleCellExperiment

#OG sce_skin command:
sce_skin <- SingleCellExperiment(assays=list(counts=counts_skin), colData=metadata_skin)

# Determine the number of clusters (nk) and the cluster names (kids), and the number of samples (ns) and sample names (sids)
colData(sce_skin)$sample_id <- as.factor(colData(sce_skin)$sample)
colData(sce_skin)$sex_id <- as.factor(colData(sce_skin)$sex)
colData(sce_skin)$cluster_id <- as.factor(colData(sce_skin)$CellType)
kids_skin <- purrr::set_names(levels(sce_skin$cluster_id))
nk_skin <- length(kids_skin)
sids_skin <- purrr::set_names(levels(sce_skin$sample_id))
ns_skin <- length(sids_skin)

# Turn named vector into a numeric vector
n_cells_skin <- as.numeric(table(sce_skin$sample_id))

# Reorder samples (rows) of the metadata to match the order of the sample names
m_skin <- match(sids_skin, sce_skin$sample_id)

# Create sample level metadata by combining the reordered metadata with the number of cells
ei_skin <- data.frame(colData(sce_skin)[m_skin, ], n_cells_skin, row.names=NULL) %>% select(-"cluster_id")

# Remove lowly expressed genes
sce_skin <- sce_skin[rowSums(counts(sce_skin)>1)>=10,]	

# Aggregate counts per sample_id and cluster_id
groups_skin <- colData(sce_skin)[, c("cluster_id", "sample_id")]
pb_skin <- aggregate.Matrix(t(counts(sce_skin)), groupings=groups_skin, fun="sum")


##########################################################
#### Turn your ase count data into a matrix:

summed_ac_counts <- cleaned_df %>%
  group_by(cells, geneid) %>%
  summarise(ref_count = sum(ref, na.rm = T), alt_count = sum(alt, na.rm = T), .groups = "drop")

any(duplicated(summed_ac_counts$cells)) #This should be TRUE

summed_ac_counts <- summed_ac_counts %>%
  mutate(cells = str_replace(cells, "^Lib19", "F1_Lib19"),
         cells = str_replace(cells, "^Lib26", "F2_Lib26"),
         cells = str_replace(cells, "^Lib34", "F3_Lib34"),
         cells = str_replace(cells, "^Lib18", "M1_Lib18"),
         cells = str_replace(cells, "^Lib32", "M2_Lib32"),
         cells = str_replace(cells, "^Lib33", "M3_Lib33"))

### Cluster your summed_ac_counts into cell-types
#All the info for the match is in groups_skins

groups_skin_name <- data.frame(groups_skin) 
groups_skin_name <- groups_skin_name %>% rownames_to_column(var = "cells")

match(groups_skin_name$cells, summed_ac_counts$cells)

group_name_ac_skin <- summed_ac_counts %>%
  left_join(groups_skin_name, by = "cells")

cleaned_group_name_ac_skin <- na.omit(group_name_ac_skin)


################ Turn into matrix ####################
############## save as dgCMatrix #####################

# Build ONE big sparse matrix (cells x genes) for ref_count & alt_count
ref_matrix_all <- cleaned_group_name_ac_skin %>%
  select(cells, geneid, ref_count) %>%
  group_by(cells, geneid) %>%
  summarise(count = sum(ref_count), .groups = "drop") %>%
  pivot_wider(names_from = geneid, values_from = count, values_fill = 0) %>%
  column_to_rownames("cells") %>%
  as.matrix() %>%
  Matrix(sparse = TRUE)


alt_matrix_all <- cleaned_group_name_ac_skin %>%
  select(cells, geneid, alt_count) %>%
  group_by(cells, geneid) %>%
  summarise(count = sum(alt_count), .groups = "drop") %>%
  pivot_wider(names_from = geneid, values_from = count, values_fill = 0) %>%
  column_to_rownames("cells") %>%
  as.matrix() %>%
  Matrix(sparse = TRUE)

cell_metadata <- cleaned_group_name_ac_skin %>%
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
ref_splitf_skin <- sapply(stringr::str_split(rownames(aggregated_ref), pattern="_", n=2), '[', 1) # Check that each cluster is present in all samples
alt_splitf_skin <- sapply(stringr::str_split(rownames(aggregated_alt), pattern="_", n=2), '[', 1) # Check that each cluster is present in all samples

ref_pb_skin <- split.data.frame(aggregated_ref, factor(ref_splitf_skin)) %>% lapply(function(u) set_colnames(t(u), stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))
alt_pb_skin <- split.data.frame(aggregated_alt, factor(alt_splitf_skin)) %>% lapply(function(u) set_colnames(t(u), stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

##############################################################################
########### Find MAR on all celltypes before splitting them up: ##############
##############################################################################

working_ac_skin <- summed_ac_counts
working_ac_skin <- working_ac_skin %>%
  filter(ref_count >= 10 & alt_count >= 10)
working_ac_skin <- working_ac_skin %>%
  mutate(sample_id = case_when( 
    grepl("^M1", cells) ~ "M1",
    grepl("^M2", cells) ~ "M2",
    grepl("^M3", cells) ~ "M3",
    grepl("^F1", cells) ~ "F1",
    grepl("^F2", cells) ~ "F2",
    grepl("^F3", cells) ~ "F3"))
working_ac_skin$log10_ref <- log10(working_ac_skin$ref_count)
working_ac_skin$log10_alt <- log10(working_ac_skin$alt_count)


#write.table(working_ac_skin, "scASE_results/skin_filtered10.txt", quote = F, row.names = F)
write.table(working_ac_skin, "cellSNP_scAlleleCount/R_results/skin_filtered10.txt", quote = F, row.names = F)

plotA <- ggplot(working_ac_skin, aes(x = sample_id, y = log10_ref)) + geom_violin(trim = F) + ggtitle("Skin") + ylab("log10 Ref Count")
plotB <- ggplot(working_ac_skin, aes(x = sample_id, y = log10_alt)) + geom_violin(trim = F) + ggtitle("Skin") + ylab("log10 Alt Count")

plotA + plotB

######################################
#### Filtering the allele counts #####
######################################

working_ac_skin %>%
  subset(sample_id == "F1") %>% pull(ref_count) %>% quantile(probs = 0.2) # 16
working_ac_skin %>%
  subset(sample_id == "F1") %>% pull(ref_count) %>% quantile(probs = 0.8) # 107
working_ac_skin %>%
  subset(sample_id == "F1") %>% pull(alt_count) %>% quantile(probs = 0.2) # 13
working_ac_skin %>%
  subset(sample_id == "F1") %>% pull(alt_count) %>% quantile(probs = 0.8) # 49 


working_ac_skin %>%
  subset(sample_id == "F2") %>% pull(ref_count) %>% quantile(probs = 0.2) # 16
working_ac_skin %>%
  subset(sample_id == "F2") %>% pull(ref_count) %>% quantile(probs = 0.8) # 116
working_ac_skin %>%
  subset(sample_id == "F2") %>% pull(alt_count) %>% quantile(probs = 0.2) # 13
working_ac_skin %>%
  subset(sample_id == "F2") %>% pull(alt_count) %>% quantile(probs = 0.8) # 52

working_ac_skin %>%
  subset(sample_id == "F3") %>% pull(ref_count) %>% quantile(probs = 0.2) # 16
working_ac_skin %>%
  subset(sample_id == "F3") %>% pull(ref_count) %>% quantile(probs = 0.8) # 110
working_ac_skin %>%
  subset(sample_id == "F3") %>% pull(alt_count) %>% quantile(probs = 0.2) # 13
working_ac_skin %>%
  subset(sample_id == "F3") %>% pull(alt_count) %>% quantile(probs = 0.8) # 51

working_ac_skin %>%
  subset(sample_id == "M1") %>% pull(ref_count) %>% quantile(probs = 0.2) # 19
working_ac_skin %>%
  subset(sample_id == "M1") %>% pull(ref_count) %>% quantile(probs = 0.8) # 155
working_ac_skin %>%
  subset(sample_id == "M1") %>% pull(alt_count) %>% quantile(probs = 0.2) # 14
working_ac_skin %>%
  subset(sample_id == "M1") %>% pull(alt_count) %>% quantile(probs = 0.8) # 62

working_ac_skin %>%
  subset(sample_id == "M2") %>% pull(ref_count) %>% quantile(probs = 0.2) # 21
working_ac_skin %>%
  subset(sample_id == "M2") %>% pull(ref_count) %>% quantile(probs = 0.8) # 204
working_ac_skin %>%
  subset(sample_id == "M2") %>% pull(alt_count) %>% quantile(probs = 0.2) # 14
working_ac_skin %>%
  subset(sample_id == "M2") %>% pull(alt_count) %>% quantile(probs = 0.8) # 75

working_ac_skin %>%
  subset(sample_id == "M3") %>% pull(ref_count) %>% quantile(probs = 0.2) # 17
working_ac_skin %>%
  subset(sample_id == "M3") %>% pull(ref_count) %>% quantile(probs = 0.8) # 112
working_ac_skin %>%
  subset(sample_id == "M3") %>% pull(alt_count) %>% quantile(probs = 0.2) # 13
working_ac_skin %>%
  subset(sample_id == "M3") %>% pull(alt_count) %>% quantile(probs = 0.8) # 52



#Pick the lowest distribution for males and highest distribution for females separately
#This is so that the allele count distribution is closest between sexes
#Less that 80% MAF for reference allele, greater than 20% MAF for alternate allele:

#ref_count > 16 female & alt_count > 13 female & alt_count < 49
#ref_count > 17 male & alt_count > 13 male & alt_count < 52

working_ac_skin_filtered <- working_ac_skin %>%
  filter(
    (sample_id == "F1" & ref_count < 53 & alt_count > 12) |
    (sample_id == "F2" & ref_count < 53 & alt_count > 12) |
    (sample_id == "F3" & ref_count < 53 & alt_count > 12) |
    (sample_id == "M1" & ref_count < 75 & alt_count > 12) |
    (sample_id == "M2" & ref_count < 75 & alt_count > 12) |
    (sample_id == "M3" & ref_count < 75 & alt_count > 12))

plotC <- ggplot(working_ac_skin_filtered, aes(x = sample_id, y = log10_ref)) + geom_violin(trim = F) + ggtitle("Skin Filtered") + ylab("log10 Ref Count")
plotD <- ggplot(working_ac_skin_filtered, aes(x = sample_id, y = log10_alt)) + geom_violin(trim = F) + ggtitle("Skin Filtered") + ylab("log10 Alt Count")

plotC + plotD

#working_ac_skin_filtered <- working_ac_skin


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

working_ac_skin_filtered <- add_chrom_column(working_ac_skin_filtered, gene_info)
working_ac_skin_filtered$chrom_type <- ifelse(working_ac_skin_filtered$chrom == "LG12", "Sex Chromosome", "Autosome")
working_ac_skin_filtered$chrom_type[working_ac_skin_filtered$geneid %in% gene_list] <- "Impacted Gene"

working_ac_skin_A <- subset(working_ac_skin_filtered, chrom_type == "Autosome")
working_ac_skin_sc <- subset(working_ac_skin_filtered, chrom_type == "Sex Chromosome")
working_ac_skin_HI <- subset(working_ac_skin_filtered, chrom_type == "Impacted Gene")

combo_working_ac_skin <- rbind(working_ac_skin_A, working_ac_skin_HI, working_ac_skin_sc)


celltype_ase <- combo_working_ac_skin %>%
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

skin_impactgenes <- subset(plot_ase, chrom_type == "Impacted Gene")
skin_autosome <- subset(plot_ase, chrom_type == "Autosome")
skin_sexchromo <- subset(plot_ase, chrom_type == "Sex Chromosome")

### Plotting:

skin_impactgenes %>% filter(sex == "F") %>% pull(avg_mar) %>% median()
skin_impactgenes %>% filter(sex == "M") %>% pull(avg_mar) %>% median()
skin_autosome %>% filter(sex == "F") %>% pull(avg_mar) %>% median()
skin_autosome %>% filter(sex == "M") %>% pull(avg_mar) %>% median()
skin_sexchromo %>% filter(sex == "F") %>% pull(avg_mar) %>% median()
skin_sexchromo %>% filter(sex == "M") %>% pull(avg_mar) %>% median()

plot_skin_IG <- ggplot(skin_impactgenes, aes(x = avg_mar)) + geom_density(aes(fill = sex), alpha = 0.5) + 
  scale_fill_manual(values = c("firebrick3", "cornflowerblue")) +
  #scale_fill_manual(values = c("firebrick3", "hotpink", "pink", "navy", "blue", "cornflowerblue")) + 
  labs(fill = "Sex") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 12)) +
  ylab("Density") + xlab("Major Allele Ratio") + ggtitle("Skin (Impacted Genes)") +
  geom_vline(aes(xintercept=0.6084237),color="firebrick4",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5651673),color="navy",linetype="dashed",size=0.8) +
  theme_classic() + theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    plot.title = element_text(size = 20, hjust = 0.5)
  )
plot_skin_IG

plot_skin_A <- ggplot(skin_autosome, aes(x = avg_mar)) + geom_density(aes(fill = sex), alpha = 0.5) + 
  scale_fill_manual(values = c("firebrick3", "cornflowerblue")) +
  #scale_fill_manual(values = c("firebrick3", "hotpink", "pink", "navy", "blue", "cornflowerblue")) + 
  labs(fill = "Sex") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 12)) +
  ylab("Density") + xlab("Major Allele Ratio") + ggtitle("Skin (Autosomes)") +
  geom_vline(aes(xintercept=0.5630952),color="firebrick4",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5709903),color="navy",linetype="dashed",size=0.8) +
  theme_classic() + theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    plot.title = element_text(size = 20, hjust = 0.5)
  )
plot_skin_A

plot_skin_SC <- ggplot(skin_sexchromo, aes(x = avg_mar)) + geom_density(aes(fill = sex), alpha = 0.5) + 
  scale_fill_manual(values = c("firebrick3", "cornflowerblue")) +
  #scale_fill_manual(values = c("firebrick3", "hotpink", "pink", "navy", "blue", "cornflowerblue")) + 
  labs(fill = "Sex") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 12)) +
  ylab("Density") + xlab("Major Allele Ratio") + ggtitle("Skin (Sex Chromosome)") +
  geom_vline(aes(xintercept=0.5663522),color="firebrick4",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5760194),color="navy",linetype="dashed",size=0.8) +
  theme_classic() + theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    plot.title = element_text(size = 20, hjust = 0.5)
  )
plot_skin_SC

skin_plots <- (plot_skin_IG | plot_skin_A | plot_skin_SC) + plot_layout(guides = "collect")
skin_plots

wilcox.test(skin_impactgenes$avg_mar ~ skin_impactgenes$sex, p.adjust.methods = "bonferroni")
wilcox.test(skin_autosome$avg_mar ~ skin_autosome$sex, p.adjust.methods = "bonferroni")
wilcox.test(skin_sexchromo$avg_mar ~ skin_sexchromo$sex, p.adjust.methods = "bonferroni")



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
write.csv(res_df_male, "scASE_results/skin_male_all.csv", row.names = T)
write.csv(sig_res_df_male, "scASE_results/skin_male_all_sig.csv", row.names = T)

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
write.csv(res_df_female, "scASE_results/skin_female_all.csv", row.names = T)
write.csv(sig_res_df_female, "scASE_results/skin_female_all_sig.csv", row.names = T)

######################## Find MAR and plot ############################

########## To directly compare celltypes ##########

# design <- ~condition + condition:sample + condition:count
# We analyzed allele-specific reads from each cell population separately with DESeq2 with the model 
# Allele-specific expression reveals genetic drivers of tissue regeneration in mice
# “~tissue + tissue:sample + tissue:allele” (where “:” denotes an interaction term between two variables in DESeq2) - from https://doi.org/10.1016/j.stem.2023.08.010
# design <- ~sex + sex:cell + sex:allele 
# This result gives me: 
# 1. baseline allele expression by sex
# 2. cell-specific effects nested within each sex
# 3. allele-specific expression interaction with sex

########### To run this for all celltypes in your matrix ##############

# Write a function so you can pull out every matrix from your ^_pb list:
# Separate out male vs. female

run_scASE_by_sex <- function(cell_type, ref_list, alt_list, output_dir = "scASE_results/skin/") {
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

cell_types <- names(ref_pb_skin)

results_by_sex <- lapply(cell_types, function(ct) {
  run_scASE_by_sex(ct, ref_pb_skin, alt_pb_skin)
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

    # 1. Epidermal

epi_F <- read.csv("scASE_results/skin/Epidermal_F.csv")
epi_M <- read.csv("scASE_results/skin/Epidermal_M.csv")
epi_F$sex = "F"
epi_M$sex = "M"
epi_combo <- rbind(epi_F, epi_M)

epi_combo$prop <- (2^epi_combo$log2FoldChange)/(1+2^epi_combo$log2FoldChange)
epi_combo$mar <- pmax(epi_combo$prop, 1 - epi_combo$prop)
epi_combo <- add_chrom_column(epi_combo, gene_info)
epi_combo$chrom_type <- ifelse(epi_combo$chrom == "LG12", "Sex Chromosome", "Autosome")
epi_combo$chrom_type[epi_combo$X %in% gene_list] <- "Impacted Gene"
epi_combo_M <- subset(epi_combo, sex == "M")
epi_combo_F <- subset(epi_combo, sex == "F")

epi_combo_M %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.6517441
epi_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.5872508
epi_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6694824


epi_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.6437316
epi_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.656608
epi_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6489218

write.csv(epi_combo_M, "scASE_results/skin/split_chromo/epidermal_M.csv", quote = FALSE, row.names = FALSE)
write.csv(epi_combo_F, "scASE_results/skin/split_chromo/epidermal_F.csv", quote = FALSE, row.names = FALSE)


plot_epi_M <- ggplot(epi_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.2) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 6)) +
  geom_vline(aes(xintercept=0.6517441),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5872508),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6694824),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Skin (Epidermal M)")
plot_epi_M

plot_epi_F <- ggplot(epi_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 2) + 
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 6)) +
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  geom_vline(aes(xintercept=0.6437316),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.656608),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6489218),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Skin (Epidermal F)")

plot_epi_F



    # 2. Fibroblast

fibro_F <- read.csv("scASE_results/skin/Fibroblast_F.csv")
fibro_M <- read.csv("scASE_results/skin/Fibroblast_M.csv")
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
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.6400939
fibro_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.6596476
fibro_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6527446

fibro_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) #  0.6144445
fibro_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.6689411
fibro_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6149924

write.csv(fibro_combo_M, "scASE_results/skin/split_chromo/fibroblast_M.csv", quote = FALSE, row.names = FALSE)
write.csv(fibro_combo_F, "scASE_results/skin/split_chromo/fibroblast_F.csv", quote = FALSE, row.names = FALSE)

plot_fibro_M <- ggplot(fibro_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 6)) +
  geom_vline(aes(xintercept=0.6400939),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6596476),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6527446),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Skin (Fibroblast M)")
plot_fibro_M

plot_fibro_F <- ggplot(fibro_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.2) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 6)) +
  geom_vline(aes(xintercept=0.6144445),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6689411),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6149924),color="darkorchid4",linetype="dashed",size=0.8) + theme_classic() + ggtitle("Skin (Fibroblast F)")
plot_fibro_F


    # 3. Granulocyte

granu_F <- read.csv("scASE_results/skin/Granulocyte_F.csv")
granu_M <- read.csv("scASE_results/skin/Granulocyte_M.csv")
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
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.6542111
granu_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.5901493
granu_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6524313

granu_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) #  0.6210538
granu_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.6925482
granu_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6466485


write.csv(granu_combo_M, "scASE_results/skin/split_chromo/granulocyte_M.csv", quote = FALSE, row.names = FALSE)
write.csv(granu_combo_F, "scASE_results/skin/split_chromo/granulocyte_F.csv", quote = FALSE, row.names = FALSE)



plot_granu_M <- ggplot(granu_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.2) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 6)) +
  geom_vline(aes(xintercept=0.6542111),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5901493),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6524313),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Skin (Granulocyte M)")
plot_granu_M 

plot_granu_F <- ggplot(granu_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.2) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 6)) +
  geom_vline(aes(xintercept=0.6210538),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6925482),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept= 0.6466485),color="darkorchid4",linetype="dashed",size=0.8) + theme_classic() + ggtitle("Skin (Granulocyte F)")
plot_granu_F



    # 4. Immune

immune_F <- read.csv("scASE_results/skin/Immune_F.csv")
immune_M <- read.csv("scASE_results/skin/Immune_M.csv")
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
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.655925
immune_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.6327101
immune_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6376234

immune_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) #  0.6413882
immune_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.5747511
immune_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6418469


write.csv(immune_combo_M, "scASE_results/skin/split_chromo/immune_M.csv", quote = FALSE, row.names = FALSE)
write.csv(immune_combo_F, "scASE_results/skin/split_chromo/immune_F.csv", quote = FALSE, row.names = FALSE)



plot_immune_M <- ggplot(immune_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.5) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 6)) +
  geom_vline(aes(xintercept=0.655925),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6327101),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept= 0.6376234),color="darkorchid4",linetype="dashed",size=0.8) +
  theme_classic() + ggtitle("Skin (Immune M)")
plot_immune_M

plot_immune_F <- ggplot(immune_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.8) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 6)) +
  geom_vline(aes(xintercept=0.6413882),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5747511),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept= 0.6418469),color="darkorchid4",linetype="dashed",size=0.8) +
  theme_classic() + ggtitle("Skin (Immune F)")
plot_immune_F


    # 5. Keratinocyte

kera_F <- read.csv("scASE_results/skin/Keratinoycte_F.csv")
kera_M <- read.csv("scASE_results/skin/Keratinoycte_M.csv")
kera_F$sex = "F"
kera_M$sex = "M"
kera_combo <- rbind(kera_F, kera_M)

kera_combo$prop <- (2^kera_combo$log2FoldChange)/(1+2^kera_combo$log2FoldChange)
kera_combo$mar <- pmax(kera_combo$prop, 1 - kera_combo$prop)
kera_combo <- add_chrom_column(kera_combo, gene_info)
kera_combo$chrom_type <- ifelse(kera_combo$chrom == "LG12", "Sex Chromosome", "Autosome")
kera_combo$chrom_type[kera_combo$X %in% gene_list] <- "Impacted Gene"
kera_combo_M <- subset(kera_combo, sex == "M")
kera_combo_F <- subset(kera_combo, sex == "F")

kera_combo_M %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.6577221
kera_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.6102204
kera_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.682996

kera_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) #  0.6184682
kera_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.6964525
kera_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6520132


write.csv(kera_combo_M, "scASE_results/skin/split_chromo/keratinocyte_M.csv", quote = FALSE, row.names = FALSE)
write.csv(kera_combo_F, "scASE_results/skin/split_chromo/keratinocyte_F.csv", quote = FALSE, row.names = FALSE)


plot_kera_M <- ggplot(kera_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.2) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 6)) +
  geom_vline(aes(xintercept=0.6577221),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6102204),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept= 0.682996),color="darkorchid4",linetype="dashed",size=0.8) +
  theme_classic() + ggtitle("Skin (Keratinocyte M)")
plot_kera_M

plot_kera_F <- ggplot(kera_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 6)) +
  geom_vline(aes(xintercept=0.6184682),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6964525),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept= 0.6520132),color="darkorchid4",linetype="dashed",size=0.8) +
  theme_classic() + ggtitle("Skin (Keratinocyte F)")
plot_kera_F


    # 6. Macrophage

macro_F <- read.csv("scASE_results/skin/Macrophage_F.csv")
macro_M <- read.csv("scASE_results/skin/Macrophage_M.csv")
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
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.6349438
macro_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.6259215
macro_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6601045

macro_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) #  0.6265234
macro_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.7079163
macro_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6754108


write.csv(macro_combo_M, "scASE_results/skin/split_chromo/macrophage_M.csv", quote = FALSE, row.names = FALSE)
write.csv(macro_combo_F, "scASE_results/skin/split_chromo/macrophage_F.csv", quote = FALSE, row.names = FALSE)



plot_macro_M <- ggplot(macro_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.2) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 6)) +
  geom_vline(aes(xintercept=0.6349438),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6259215),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept= 0.6601045),color="darkorchid4",linetype="dashed",size=0.8) +
  theme_classic() + ggtitle("Skin (Macrophage M)")
plot_macro_M

plot_macro_F <- ggplot(macro_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.2) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 6)) +
  geom_vline(aes(xintercept=0.6265234),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.7079163),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept= 0.6754108),color="darkorchid4",linetype="dashed",size=0.8) +
  theme_classic() + ggtitle("Skin (Macrophage F)")
plot_macro_F


    # 7. Melanocyte

mela_F <- read.csv("scASE_results/skin/Melanocyte_F.csv")
mela_M <- read.csv("scASE_results/skin/Melanocyte_M.csv")
mela_F$sex = "F"
mela_M$sex = "M"
mela_combo <- rbind(mela_F, mela_M)

mela_combo$prop <- (2^mela_combo$log2FoldChange)/(1+2^mela_combo$log2FoldChange)
mela_combo$mar <- pmax(mela_combo$prop, 1 - mela_combo$prop)
mela_combo <- add_chrom_column(mela_combo, gene_info)
mela_combo$chrom_type <- ifelse(mela_combo$chrom == "LG12", "Sex Chromosome", "Autosome")
mela_combo$chrom_type[mela_combo$X %in% gene_list] <- "Impacted Gene"
mela_combo_M <- subset(mela_combo, sex == "M")
mela_combo_F <- subset(mela_combo, sex == "F")

mela_combo_M %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.6399816
mela_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.7924164
mela_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6765502

mela_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) #  0.6526242
mela_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.6984933
mela_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6737546



write.csv(mela_combo_M, "scASE_results/skin/split_chromo/melanocyte_M.csv", quote = FALSE, row.names = FALSE)
write.csv(mela_combo_F, "scASE_results/skin/split_chromo/melanocyte_F.csv", quote = FALSE, row.names = FALSE)


plot_mela_M <- ggplot(mela_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.8) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 6)) +
  geom_vline(aes(xintercept=0.6399816),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.7924164),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept= 0.6765502),color="darkorchid4",linetype="dashed",size=0.8) +
  theme_classic() + ggtitle("Skin (Melanocyte M)")
plot_mela_M

plot_mela_F <- ggplot(mela_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 6)) +
  geom_vline(aes(xintercept=0.6526242),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6984933),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept= 0.6737546),color="darkorchid4",linetype="dashed",size=0.8) +
  theme_classic() + ggtitle("Skin (Melanocyte F)")
plot_mela_F


    # 8. Mesenchymal stromal

messtro_F <- read.csv("scASE_results/skin/Mesenchymal stromal_F.csv")
messtro_M <- read.csv("scASE_results/skin/Mesenchymal stromal_M.csv")
messtro_F$sex = "F"
messtro_M$sex = "M"
messtro_combo <- rbind(messtro_F, messtro_M)

messtro_combo$prop <- (2^messtro_combo$log2FoldChange)/(1+2^messtro_combo$log2FoldChange)
messtro_combo$mar <- pmax(messtro_combo$prop, 1 - messtro_combo$prop)
messtro_combo <- add_chrom_column(messtro_combo, gene_info)
messtro_combo$chrom_type <- ifelse(messtro_combo$chrom == "LG12", "Sex Chromosome", "Autosome")
messtro_combo$chrom_type[messtro_combo$X %in% gene_list] <- "Impacted Gene"
messtro_combo_M <- subset(messtro_combo, sex == "M")
messtro_combo_F <- subset(messtro_combo, sex == "F")


messtro_combo_M %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.6457923
messtro_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.6411432
messtro_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6594207

messtro_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) #  0.6040049
messtro_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.5611245
messtro_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6010586



write.csv(messtro_combo_M, "scASE_results/skin/split_chromo/messenchymal_stromal_M.csv", quote = FALSE, row.names = FALSE)
write.csv(messtro_combo_F, "scASE_results/skin/split_chromo/messenchymal_stromal_F.csv", quote = FALSE, row.names = FALSE)



plot_messtro_M <- ggplot(messtro_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 6)) +
  geom_vline(aes(xintercept=0.6457923),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6411432),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept= 0.6594207),color="darkorchid4",linetype="dashed",size=0.8) +
  theme_classic() + ggtitle("Skin (Mesenchymal stromal M)")
plot_messtro_M

plot_messtro_F <- ggplot(messtro_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.5) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 6)) +
  geom_vline(aes(xintercept=0.6040049),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5611245),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept= 0.6010586),color="darkorchid4",linetype="dashed",size=0.8) +
  theme_classic() + ggtitle("Skin (Mesenchymal stromal F)")
plot_messtro_F


    # 9. Mitochondrial Rich

mito_F <- read.csv("scASE_results/skin/Mitochondrial rich_F.csv")
mito_M <- read.csv("scASE_results/skin/Mitochondrial rich_M.csv")
mito_F$sex = "F"
mito_M$sex = "M"
mito_combo <- rbind(mito_F, mito_M)

mito_combo$prop <- (2^mito_combo$log2FoldChange)/(1+2^mito_combo$log2FoldChange)
mito_combo$mar <- pmax(mito_combo$prop, 1 - mito_combo$prop)
mito_combo <- add_chrom_column(mito_combo, gene_info)
mito_combo$chrom_type <- ifelse(mito_combo$chrom == "LG12", "Sex Chromosome", "Autosome")
mito_combo$chrom_type[mito_combo$X %in% gene_list] <- "Impacted Gene"
mito_combo_M <- subset(mito_combo, sex == "M")
mito_combo_F <- subset(mito_combo, sex == "F")


mito_combo_M %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.6596609
mito_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.7221019
mito_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.7068482

mito_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) #  0.6265359
mito_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.6109287
mito_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6402949


write.csv(mito_combo_M, "scASE_results/skin/split_chromo/mitochondrial_M.csv", quote = FALSE, row.names = FALSE)
write.csv(mito_combo_F, "scASE_results/skin/split_chromo/mitochondrial_F.csv", quote = FALSE, row.names = FALSE)

plot_mito_M <- ggplot(mito_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 6)) +
  geom_vline(aes(xintercept=0.6596609),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.7221019),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept= 0.7068482),color="darkorchid4",linetype="dashed",size=0.8) +
  theme_classic() + ggtitle("Skin (Mitochondrial rich M)")
plot_mito_M

plot_mito_F <- ggplot(mito_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.5) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 6)) +
  geom_vline(aes(xintercept=0.6265359),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6109287),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept= 0.6402949),color="darkorchid4",linetype="dashed",size=0.8) +
  theme_classic() + ggtitle("Skin (Mitochondrial rich F)")
plot_mito_F


    # 10. Peridermal

peri_F <- read.csv("scASE_results/skin/Peridermal_F.csv")
peri_M <- read.csv("scASE_results/skin/Peridermal_M.csv")
peri_F$sex = "F"
peri_M$sex = "M"
peri_combo <- rbind(peri_F, peri_M)

peri_combo$prop <- (2^peri_combo$log2FoldChange)/(1+2^peri_combo$log2FoldChange)
peri_combo$mar <- pmax(peri_combo$prop, 1 - peri_combo$prop)
peri_combo <- add_chrom_column(peri_combo, gene_info)
peri_combo$chrom_type <- ifelse(peri_combo$chrom == "LG12", "Sex Chromosome", "Autosome")
peri_combo$chrom_type[peri_combo$X %in% gene_list] <- "Impacted Gene"
peri_combo_M <- subset(peri_combo, sex == "M")
peri_combo_F <- subset(peri_combo, sex == "F")


peri_combo_M %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.6337985
peri_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.6601617
peri_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6426094

peri_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) #  0.6101088
peri_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) #  0.6510848
peri_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6668388

write.csv(peri_combo_M, "scASE_results/skin/split_chromo/peridermal_M.csv", quote = FALSE, row.names = FALSE)
write.csv(peri_combo_F, "scASE_results/skin/split_chromo/peridermal_F.csv", quote = FALSE, row.names = FALSE)


plot_peri_M <- ggplot(peri_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 6)) +
  geom_vline(aes(xintercept=0.6337985),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6601617),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept= 0.6426094),color="darkorchid4",linetype="dashed",size=0.8) +
  theme_classic() + ggtitle("Skin (Peridermal M)")
plot_peri_M

plot_peri_F <- ggplot(peri_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.5) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 6)) +
  geom_vline(aes(xintercept=0.6101088),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6510848),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept= 0.6668388),color="darkorchid4",linetype="dashed",size=0.8) +
  theme_classic() + ggtitle("Skin (Peridermal F)")
plot_peri_F


    # 11. Pigment

pig_F <- read.csv("scASE_results/skin/Pigment_F.csv")
pig_M <- read.csv("scASE_results/skin/Pigment_M.csv")
pig_F$sex = "F"
pig_M$sex = "M"
pig_combo <- rbind(pig_F, pig_M)

pig_combo$prop <- (2^pig_combo$log2FoldChange)/(1+2^pig_combo$log2FoldChange)
pig_combo$mar <- pmax(pig_combo$prop, 1 - pig_combo$prop)
pig_combo <- add_chrom_column(pig_combo, gene_info)
pig_combo$chrom_type <- ifelse(pig_combo$chrom == "LG12", "Sex Chromosome", "Autosome")
pig_combo$chrom_type[pig_combo$X %in% gene_list] <- "Impacted Gene"
pig_combo_M <- subset(pig_combo, sex == "M")
pig_combo_F <- subset(pig_combo, sex == "F")

pig_combo_M %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.6609114
pig_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.6735448
pig_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.7269015

pig_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.6184696
pig_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.6694827
pig_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6287761

write.csv(pig_combo_M, "scASE_results/skin/split_chromo/pigment_M.csv", quote = FALSE, row.names = FALSE)
write.csv(pig_combo_F, "scASE_results/skin/split_chromo/pigment_F.csv", quote = FALSE, row.names = FALSE)

plot_pig_M <- ggplot(pig_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 3) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 6)) +
  geom_vline(aes(xintercept=0.6609114),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6735448),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept= 0.7269015),color="darkorchid4",linetype="dashed",size=0.8) +
  theme_classic() + ggtitle("Skin (Pigment M)")
plot_pig_M

plot_pig_F <- ggplot(pig_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 2.4) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 6)) +
  geom_vline(aes(xintercept=0.6184696),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6694827),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept= 0.6287761),color="darkorchid4",linetype="dashed",size=0.8) +
  theme_classic() + ggtitle("Skin (Pigment F)")
plot_pig_F



    # 12. Ribosomal Protein Rich

ribo_F <- read.csv("scASE_results/skin/Ribosomal protein rich_F.csv")
ribo_M <- read.csv("scASE_results/skin/Ribosomal protein rich_M.csv")
ribo_F$sex = "F"
ribo_M$sex = "M"
ribo_combo <- rbind(ribo_F, ribo_M)

ribo_combo$prop <- (2^ribo_combo$log2FoldChange)/(1+2^ribo_combo$log2FoldChange)
ribo_combo$mar <- pmax(ribo_combo$prop, 1 - ribo_combo$prop)
ribo_combo <- add_chrom_column(ribo_combo, gene_info)
ribo_combo$chrom_type <- ifelse(ribo_combo$chrom == "LG12", "Sex Chromosome", "Autosome")
ribo_combo$chrom_type[ribo_combo$X %in% gene_list] <- "Impacted Gene"
ribo_combo_M <- subset(ribo_combo, sex == "M")
ribo_combo_F <- subset(ribo_combo, sex == "F")

ribo_combo_M %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.6189474
ribo_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.6073443
ribo_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6270835

ribo_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.6023172
ribo_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.647777
ribo_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6208195

write.csv(ribo_combo_M, "scASE_results/skin/split_chromo/ribosomal_M.csv", quote = FALSE, row.names = FALSE)
write.csv(ribo_combo_F, "scASE_results/skin/split_chromo/ribosomal_F.csv", quote = FALSE, row.names = FALSE)

plot_ribo_M <- ggplot(ribo_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 2) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 6)) +
  geom_vline(aes(xintercept=0.6189474),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6073443),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept= 0.6270835),color="darkorchid4",linetype="dashed",size=0.8) +
  theme_classic() + ggtitle("Skin (Ribosomal protein rich M)")
plot_ribo_M

plot_ribo_F <- ggplot(ribo_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 2.5) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 6)) +
  geom_vline(aes(xintercept=0.6023172),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.647777),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept= 0.6208195),color="darkorchid4",linetype="dashed",size=0.8) +
  theme_classic() + ggtitle("Skin (Ribosomal protein rich F)")
plot_ribo_F



    # 13. Skeletal Muscle

skel_F <- read.csv("scASE_results/skin/Skeletal muscle_F.csv")
skel_M <- read.csv("scASE_results/skin/Skeletal muscle_M.csv")
skel_F$sex = "F"
skel_M$sex = "M"
skel_combo <- rbind(skel_F, skel_M)

skel_combo$prop <- (2^skel_combo$log2FoldChange)/(1+2^skel_combo$log2FoldChange)
skel_combo$mar <- pmax(skel_combo$prop, 1 - skel_combo$prop)
skel_combo <- add_chrom_column(skel_combo, gene_info)
skel_combo$chrom_type <- ifelse(skel_combo$chrom == "LG12", "Sex Chromosome", "Autosome")
skel_combo$chrom_type[skel_combo$X %in% gene_list] <- "Impacted Gene"
skel_combo_M <- subset(skel_combo, sex == "M")
skel_combo_F <- subset(skel_combo, sex == "F")


skel_combo_M %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.7187436
skel_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.5911652
skel_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.7324035

skel_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.7175138
skel_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.6846121
skel_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.7248725

write.csv(skel_combo_M, "scASE_results/skin/split_chromo/skeletal_M.csv", quote = FALSE, row.names = FALSE)
write.csv(skel_combo_F, "scASE_results/skin/split_chromo/skeletal_F.csv", quote = FALSE, row.names = FALSE)


plot_skel_M <- ggplot(skel_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.2) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 6)) +
  geom_vline(aes(xintercept=0.7187436),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5911652),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept= 0.7324035),color="darkorchid4",linetype="dashed",size=0.8) +
  theme_classic() + ggtitle("Skin (Skeletal muscle M)")
plot_skel_M

plot_skel_F <- ggplot(skel_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.2) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 6)) +
  geom_vline(aes(xintercept=0.7175138),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6846121),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept= 0.7248725),color="darkorchid4",linetype="dashed",size=0.8) +
  theme_classic() + ggtitle("Skin (Skeletal muscle)")
plot_skel_F

    # 14. Stromal

strom_F <- read.csv("scASE_results/skin/Stromal_F.csv")
strom_M <- read.csv("scASE_results/skin/Stromal_M.csv")
strom_F$sex = "F"
strom_M$sex = "M"
strom_combo <- rbind(strom_F, strom_M)

strom_combo$prop <- (2^strom_combo$log2FoldChange)/(1+2^strom_combo$log2FoldChange)
strom_combo$mar <- pmax(strom_combo$prop, 1 - strom_combo$prop)
strom_combo <- add_chrom_column(strom_combo, gene_info)
strom_combo$chrom_type <- ifelse(strom_combo$chrom == "LG12", "Sex Chromosome", "Autosome")
strom_combo$chrom_type[strom_combo$X %in% gene_list] <- "Impacted Gene"
strom_combo_M <- subset(strom_combo, sex == "M")
strom_combo_F <- subset(strom_combo, sex == "F")

strom_combo_M %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) #  0.6299686
strom_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.6110444
strom_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6302789

strom_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.6159705
strom_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.6406814
strom_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6097933

write.csv(strom_combo_M, "scASE_results/skin/split_chromo/stromal_M.csv", quote = FALSE, row.names = FALSE)
write.csv(strom_combo_F, "scASE_results/skin/split_chromo/stromal_F.csv", quote = FALSE, row.names = FALSE)


plot_strom_M <- ggplot(strom_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.5) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 6)) +
  geom_vline(aes(xintercept=0.6299686),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6110444),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept= 0.6302789),color="darkorchid4",linetype="dashed",size=0.8) +
  theme_classic() + ggtitle("Skin (Stromal M)")
plot_strom_M

plot_strom_F <- ggplot(strom_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.5) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 6)) +
  geom_vline(aes(xintercept=0.6159705),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6406814),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept= 0.6097933),color="darkorchid4",linetype="dashed",size=0.8) +
  theme_classic() + ggtitle("Skin (Stromal F)")
plot_strom_F


#### Combo your plots ####

library(patchwork)

skin_celltype_plots_M <- wrap_plots(list(plot_epi_M, plot_fibro_M, plot_granu_M, plot_immune_M, plot_kera_M, 
                                       plot_macro_M, plot_mela_M, plot_messtro_M, plot_mito_M, plot_peri_M,
                                       plot_pig_M, plot_ribo_M, plot_skel_M, plot_strom_M), ncol = 3, nrow = 5) +
  plot_layout(guides = "collect")

skin_celltype_plots_M

skin_celltype_plots_F <- wrap_plots(list(plot_epi_F, plot_fibro_F, plot_granu_F, plot_immune_F, plot_kera_F, 
                                         plot_macro_F, plot_mela_F, plot_messtro_F, plot_mito_F, plot_peri_F,
                                         plot_pig_F, plot_ribo_F, plot_skel_F, plot_strom_F), ncol = 3, nrow = 5) +
  plot_layout(guides = "collect")

skin_celltype_plots_F



######################################################
######## Combine celltypes to test for ASE ###########
######################################################


all_celltype_M <- rbind(epi_combo_M, fibro_combo_M, granu_combo_M, immune_combo_M, kera_combo_M, macro_combo_M,
                        mela_combo_M, messtro_combo_M, mito_combo_M, peri_combo_M, pig_combo_M, ribo_combo_M, 
                        skel_combo_M, strom_combo_M)

all_celltype_F <- rbind(epi_combo_F, fibro_combo_F, granu_combo_F, immune_combo_F, kera_combo_F, macro_combo_F,
                        mela_combo_F, messtro_combo_F, mito_combo_F, peri_combo_F, pig_combo_F, ribo_combo_F, 
                        skel_combo_F, strom_combo_F)


all_celltype_M_means <- all_celltype_M %>%
  group_by(X, chrom_type) %>%
  summarise(mar_mean = mean(mar), .groups = "drop")

all_celltype_F_means <- all_celltype_F %>%
  group_by(X, chrom_type) %>%
  summarise(mar_mean = mean(mar), .groups = "drop")



all_celltype_M_means %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar_mean)) #  0.658
all_celltype_M_means %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar_mean)) # 0.645
all_celltype_M_means %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar_mean)) # 0.668

all_celltype_F_means %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar_mean)) # 0.631
all_celltype_F_means %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar_mean)) # 0.649
all_celltype_F_means %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar_mean)) # 0.642



plot_celltype_M <- ggplot(all_celltype_M_means, aes(x = mar_mean)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.2) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 15)) +
  geom_vline(aes(xintercept=0.658),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.645),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept= 0.668),color="darkorchid4",linetype="dashed",size=0.8) +
  theme_classic() + ggtitle("Skin M (all celltypes)")
plot_celltype_M



plot_celltype_F <- ggplot(all_celltype_F_means, aes(x = mar_mean)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.2) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 15)) +
  geom_vline(aes(xintercept=0.631),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.649),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept= 0.642),color="darkorchid4",linetype="dashed",size=0.8) +
  theme_classic() + ggtitle("Skin F (all celltypes)")
plot_celltype_F


combo_everything <- wrap_plots(plot_celltype_M, plot_celltype_F) +
  plot_layout(guides = "collect")

combo_everything


###################################################
######## OLD (test for sex interaction) ###########
###################################################

#To include the allele counts:

sce_skin_AC <- SingleCellExperiment(
  assays = list(ref_counts = ref_mat, alt_counts = alt_mat),
  colData = metadata_skin
)

#################################################
################# Epidermal #####################
#################################################

epidermal_mat <- ref_pb_skin[["Epidermal"]]
dense_epidermal <- as.matrix(epidermal_mat)
epidermal_df <- as.data.frame(dense_epidermal)
epidermal_df <- tibble::rownames_to_column(epidermal_df, var = "gene_id")
epidermal_df_long <- epidermal_df %>%
  tidyr::pivot_longer(-gene_id, names_to = "cell", values_to = "ref")

epidermal_alt <- alt_pb_skin[["Epidermal"]]
dense_epi_alt <- as.matrix(epidermal_alt)
epidermal_alt_df <- as.data.frame(dense_epi_alt)
epidermal_alt_df <- tibble::rownames_to_column(epidermal_alt_df, var = "gene_id")
epidermal_alt_long <- epidermal_alt_df %>%
  tidyr::pivot_longer(-gene_id, names_to = "cell", values_to = "alt")

epidermal_df_long$alt <- epidermal_alt_long$alt

cleaned_epidermal_df_long <- epidermal_df_long %>%
  pivot_longer(cols = c(ref, alt), names_to = "allele", values_to = "count") %>%
  mutate(allele = factor(allele, levels = c("ref", "alt")))

coldata <- cleaned_epidermal_df_long %>%
  distinct(cell, allele) %>%
  mutate(cell = factor(cell), allele = factor(allele))

count_mat <- cleaned_epidermal_df_long %>%
  unite(sample, cell, allele, sep = "_") %>%
  pivot_wider(names_from = sample, values_from = count, values_fill = 0) %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

coldata <- data.frame(
  sample = colnames(count_mat),
  cell = sub("_(ref|alt)", "", colnames(count_mat)),
  allele = sub(".*(ref|alt)", "\\1", colnames(count_mat))
) #coldata rows match count matrix columns

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = coldata,
  design = ~ cell + allele
)

dds <- DESeq(dds)
resultsNames(dds)


epidermal_res <- results(dds, contrast = c("allele", "ref", "alt"))
# Genes with adjusted p-value < 0.05 and log2FoldChange ≠ 0 show evidence of allelic imbalance.
# Positive log2FoldChange means ref > alt; negative means alt > ref.
int_epidermal_res <- as.data.frame(epidermal_res)
sig_res_epidermal <- int_epidermal_res %>%
  filter(padj < 0.05, log2FoldChange != 0) 

write.csv(res_epidermal, "scASE_results/epidermal.csv", quote = FALSE, row.names = FALSE)



#################################################
################# Fibroblast ####################
#################################################

Fibroblast_mat <- ref_pb_skin[["Fibroblast"]]
dense_Fibroblast <- as.matrix(Fibroblast_mat)
Fibroblast_df <- as.data.frame(dense_Fibroblast)
Fibroblast_df <- tibble::rownames_to_column(Fibroblast_df, var = "gene_id")
Fibroblast_df_long <- Fibroblast_df %>%
  tidyr::pivot_longer(-gene_id, names_to = "cell", values_to = "ref")

Fibroblast_alt <- alt_pb_skin[["Fibroblast"]]
dense_epi_alt <- as.matrix(Fibroblast_alt)
Fibroblast_alt_df <- as.data.frame(dense_epi_alt)
Fibroblast_alt_df <- tibble::rownames_to_column(Fibroblast_alt_df, var = "gene_id")
Fibroblast_alt_long <- Fibroblast_alt_df %>%
  tidyr::pivot_longer(-gene_id, names_to = "cell", values_to = "alt")

Fibroblast_df_long$alt <- Fibroblast_alt_long$alt
Fibroblast_df_long$sex <- ifelse(grepl("F", Fibroblast_df_long$cell), "F", "M")

cleaned_Fibroblast_df_long <- Fibroblast_df_long %>%
  pivot_longer(cols = c(ref, alt), names_to = "allele", values_to = "count") %>%
  mutate(allele = factor(allele, levels = c("ref", "alt")))

cleaned_Fibroblast_df_long <- cleaned_Fibroblast_df_long %>%
  mutate(cell = factor(cell), allele = factor(allele, levels = c("ref", "alt")),
         sex = factor(sex))

cleaned_Fibroblast_df_long <- cleaned_Fibroblast_df_long %>%
  unite(sample, cell, allele, sep = "_", remove = F)

nona_Fibroblast_df_long <- na.omit(cleaned_Fibroblast_df_long)

count_mat <- cleaned_Fibroblast_df_long %>%
  select(gene_id, sample, count) %>%
  pivot_wider(names_from = sample, values_from = count, values_fill = 0) %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

coldata <- cleaned_Fibroblast_df_long %>%
  distinct(sample, cell, allele, sex) %>%
  as.data.frame()
rownames(coldata) <- coldata$sample

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = coldata,
  #design = ~ sex + sex:cell + sex:allele
  design = ~ allele
)

dds <- DESeq(dds)

#To view names:
resultsNames(dds)

Fibroblast_res <- results(dds, name = "allele_alt_vs_ref")
# Genes with adjusted p-value < 0.05 and log2FoldChange ≠ 0 show evidence of allelic imbalance.
# Positive log2FC: Alt allele is more expressed in males than in females.
# Negative log2FC: Alt allele is more expressed in females than in males

fibro.alt.vs.ref <- 2^Fibroblast_res$log2FoldChange
plot(fibro.alt.vs.ref)

int_Fibroblast_res <- as.data.frame(Fibroblast_res)
res_Fibroblast <- int_Fibroblast_res %>%
  filter(log2FoldChange != 0) 
sig_res_Fibroblast <- int_Fibroblast_res %>%
  filter(padj < 0.05, log2FoldChange != 0) 

write.csv(res_Fibroblast, "scASE_results/Fibroblast.csv", quote = FALSE, row.names = FALSE)

#################################################

# Write a function so you can pull out every matrix from your ^_pb list:

run_scASE_analysis <- function(cell_type, ref_list, alt_list, output_dir = "scASE_results/") {
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
  
  # Merge ref and alt counts
  combined_df <- ref_df %>%
    mutate(alt = alt_df$alt,
           sex = ifelse(grepl("F", cell), "F", "M")) # Adjust this line if needed
  
  # Reshape into long count table
  cleaned_long <- combined_df %>%
    pivot_longer(cols = c(ref, alt), names_to = "allele", values_to = "count") %>%
    mutate(
      cell = factor(cell),
      allele = factor(allele, levels = c("ref", "alt")),
      sex = factor(sex)
    ) %>%
    unite(sample, cell, allele, sep = "_", remove = FALSE) %>%
    na.omit()
  
  # Build count matrix
  count_mat <- cleaned_long %>%
    select(gene_id, sample, count) %>%
    pivot_wider(names_from = sample, values_from = count, values_fill = 0) %>%
    column_to_rownames("gene_id") %>%
    as.matrix()
  
  # Build colData
  coldata <- cleaned_long %>%
    distinct(sample, cell, allele, sex) %>%
    as.data.frame()
  rownames(coldata) <- coldata$sample
  
  # DESeq2
  dds <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData = coldata,
    design = ~ sex + allele + sex:allele
  )
  
  dds <- DESeq(dds)
  res <- results(dds, name = "sexM.allelealt")
  res_df <- as.data.frame(res) %>%
    filter(log2FoldChange != 0)
  
  # Output
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  write.csv(res_df, file.path(output_dir, paste0(cell_type, ".csv")), row.names = FALSE)
  
  return(res_df)
}


cell_types <- names(ref_pb_skin)
results_list <- lapply(cell_types, function(ct) {
  run_scASE_analysis(cell_type = ct, ref_list = ref_pb_skin, alt_list = alt_pb_skin)
})
names(results_list) <- cell_types

