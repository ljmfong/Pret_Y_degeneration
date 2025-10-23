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

remove_dense_snps_fast <- function(df, window_bp = 10000, snp_threshold = 500) {
  setDT(df)
  df <- df[order(CHROM, geneid, POS)]
  
  df[, snps_in_window := {
    # For each gene, use a moving window to count SNPs efficiently
    start <- 1L
    count <- integer(.N)
    for (i in seq_len(.N)) {
      while (POS[i] - POS[start] > window_bp) start <- start + 1L
      count[i] <- i - start + 1L
    }
    count
  }, by = .(CHROM, geneid)]
  
  df[snps_in_window <= snp_threshold]
}


filt_cleaned_df <- remove_dense_snps_fast(cleaned_df)


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

summed_ac_counts <- filt_cleaned_df %>%
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

write.table(working_ac_liver, "cellSNP_scAlleleCount/R_results/snpfilt_liver_filtered10.txt", quote = F, row.names = F)
