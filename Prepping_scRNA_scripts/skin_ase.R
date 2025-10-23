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
library(data.table)
library(DESeq2)
library(dplyr)
library(ggbreak)
library(ggplot2)
library(magrittr)
library(Matrix)
library(Matrix.utils)
library(purrr)
library(RcpRoll)
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

write.table(working_ac_skin, "cellSNP_scAlleleCount/R_results/skin_filtered10.txt", quote = F, row.names = F)


