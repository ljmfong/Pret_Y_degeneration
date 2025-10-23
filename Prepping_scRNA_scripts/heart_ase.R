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
library(data.table)
library(DESeq2)
library(dplyr)
library(ggbreak)
library(ggplot2)
library(magrittr)
library(Matrix)
library(Matrix.utils)
library(purrr)
library(RCurl)
library(RcpRoll)
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
write.table(working_ac_heart, "cellSNP_scAlleleCount/R_results/heart_filtered10.txt", quote = F, row.names = F)

                                                                                      
