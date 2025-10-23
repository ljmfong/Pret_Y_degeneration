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
library(data.table)
library(dplyr)
library(ggbreak)
library(ggplot2)
library(magrittr)
library(Matrix)
library(Matrix.utils)
library(purrr)
library(RCurl)
library(RcppRoll)
library(scales)
library(scuttle)
library(Seurat)
library(SingleCellExperiment)
library(stringr)
library(tibble)
library(tidyverse)

###################################################

##### 2. Load in Iulia's matrices ####

load("iulia_scRNA/gonad_umap_3_rename_clusters.RData")

counts_gonad <- gonad_umap_3_rename_clusters@assays$RNA@counts
metadata_gonad <- gonad_umap_3_rename_clusters@meta.data
metadata_gonad$CellType <- factor(gonad_umap_3_rename_clusters@active.ident)

#allelecount_Lib11 <- read.table("scrna_scAlleleCount/complete_matrice_for_R/Lib11_counts.txt", sep = " ", header = T)
#allelecount_Lib23 <- read.table("scrna_scAlleleCount/complete_matrice_for_R/Lib23_counts.txt", sep = " ", header = T)
#allelecount_Lib25 <- read.table("scrna_scAlleleCount/complete_matrice_for_R/Lib25_counts.txt", sep = " ", header = T)
#allelecount_Lib35 <- read.table("scrna_scAlleleCount/complete_matrice_for_R/Lib35_counts.txt", sep = " ", header = T)
#allelecount_Lib36 <- read.table("scrna_scAlleleCount/complete_matrice_for_R/Lib36_counts.txt", sep = " ", header = T)
#allelecount_Lib37 <- read.table("scrna_scAlleleCount/complete_matrice_for_R/Lib37_counts.txt", sep = " ", header = T)

allelecount_Lib11 <- read.table("cellSNP_scAlleleCount/complete_matrice_for_R/Lib11_counts.txt", sep = " ", header = T)
allelecount_Lib23 <- read.table("cellSNP_scAlleleCount/complete_matrice_for_R/Lib23_counts.txt", sep = " ", header = T)
allelecount_Lib25 <- read.table("cellSNP_scAlleleCount/complete_matrice_for_R/Lib25_counts.txt", sep = " ", header = T)
allelecount_Lib35 <- read.table("cellSNP_scAlleleCount/complete_matrice_for_R/Lib35_counts.txt", sep = " ", header = T)
allelecount_Lib36 <- read.table("cellSNP_scAlleleCount/complete_matrice_for_R/Lib36_counts.txt", sep = " ", header = T)
allelecount_Lib37 <- read.table("cellSNP_scAlleleCount/complete_matrice_for_R/Lib37_counts.txt", sep = " ", header = T)

merged_gonad_AC <- rbind(allelecount_Lib11, allelecount_Lib23, allelecount_Lib25, 
                          allelecount_Lib35, allelecount_Lib36, allelecount_Lib37)
colnames(merged_gonad_AC)[3] <- "cells"


#### Attach gene info and find matching gene_name for positional information:

gene_info <- read.table("new_Ret_SNPs/gene_boundary_w_gene_names.bed", header = T, sep = "\t")

attached_gene_AC <- merged_gonad_AC %>%
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
#write.table(cleaned_df, "scASE_results/gonad_cleaned_df_gene_attached.txt", quote = F, row.names = F)
write.table(cleaned_df, "cellSNP_scAlleleCount/R_results/gonad_cleaned_df_gene_attached.txt", quote = F, row.names = F)
#If loading in:
#cleaned_df <- read.table(file = "scASE_results/gonad_cleaned_df_gene_attached.txt", header = T)
cleaned_df <- read.table(file = "cellSNP_scAlleleCount/R_results/gonad_cleaned_df_gene_attached.txt", header = T)

######### Attach to your SingleCellExperiment
#OG sce_gonad command:
sce_gonad <- SingleCellExperiment(assays=list(counts=counts_gonad), colData=metadata_gonad)

# Determine the number of clusters (nk) and the cluster names (kids), and the number of samples (ns) and sample names (sids)
colData(sce_gonad)$sample_id <- as.factor(colData(sce_gonad)$sample)
colData(sce_gonad)$sex_id <- as.factor(colData(sce_gonad)$sex)
colData(sce_gonad)$cluster_id <- as.factor(colData(sce_gonad)$CellType)
kids_gonad <- purrr::set_names(levels(sce_gonad$cluster_id))
nk_gonad <- length(kids_gonad)
sids_gonad <- purrr::set_names(levels(sce_gonad$sample_id))
ns_gonad <- length(sids_gonad)

# Turn named vector into a numeric vector
n_cells_gonad <- as.numeric(table(sce_gonad$sample_id))

# Reorder samples (rows) of the metadata to match the order of the sample names
m_gonad <- match(sids_gonad, sce_gonad$sample_id)

# Create sample level metadata by combining the reordered metadata with the number of cells
ei_gonad <- data.frame(colData(sce_gonad)[m_gonad, ], n_cells_gonad, row.names=NULL) %>% select(-"cluster_id")

# Remove lowly expressed genes
sce_gonad <- sce_gonad[rowSums(counts(sce_gonad)>1)>=10,]	

# Aggregate counts per sample_id and cluster_id
groups_gonad <- colData(sce_gonad)[, c("cluster_id", "sample_id")]
pb_gonad <- aggregate.Matrix(t(counts(sce_gonad)), groupings=groups_gonad, fun="sum")


##############

summed_ac_counts <- cleaned_df %>%
  group_by(cells, geneid) %>%
  summarise(ref_count = sum(ref, na.rm = T), alt_count = sum(alt, na.rm = T), .groups = "drop")

any(duplicated(summed_ac_counts$cells)) #This should be TRUE


summed_ac_counts <- summed_ac_counts %>%
  mutate(cells = str_replace(cells, "^Lib11", "F1_Lib11"),
         cells = str_replace(cells, "^Lib23", "F2_Lib23"),
         cells = str_replace(cells, "^Lib25", "F3_Lib25"),
         cells = str_replace(cells, "^Lib35", "M4_Lib35"),
         cells = str_replace(cells, "^Lib36", "M5_Lib36"),
         cells = str_replace(cells, "^Lib37", "M6_Lib37"))

### Cluster your summed_ac_counts into cell-types
#All the info for the match is in groups_gonads

groups_gonad_name <- data.frame(groups_gonad) 
groups_gonad_name <- groups_gonad_name %>% rownames_to_column(var = "cells")

match(groups_gonad_name$cells, summed_ac_counts$cells)

group_name_ac_gonad <- summed_ac_counts %>%
  left_join(groups_gonad_name, by = "cells")

cleaned_group_name_ac_gonad <- na.omit(group_name_ac_gonad)

################ Turn into matrix ####################
############## save as dgCMatrix #####################
# Build ONE big sparse matrix (cells x genes) for ref_count & alt_count
ref_matrix_all <- cleaned_group_name_ac_gonad %>%
  select(cells, geneid, ref_count) %>%
  group_by(cells, geneid) %>%
  summarise(count = sum(ref_count), .groups = "drop") %>%
  pivot_wider(names_from = geneid, values_from = count, values_fill = 0) %>%
  column_to_rownames("cells") %>%
  as.matrix() %>%
  Matrix(sparse = TRUE)


alt_matrix_all <- cleaned_group_name_ac_gonad %>%
  select(cells, geneid, alt_count) %>%
  group_by(cells, geneid) %>%
  summarise(count = sum(alt_count), .groups = "drop") %>%
  pivot_wider(names_from = geneid, values_from = count, values_fill = 0) %>%
  column_to_rownames("cells") %>%
  as.matrix() %>%
  Matrix(sparse = TRUE)

cell_metadata <- cleaned_group_name_ac_gonad %>%
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

ref_splitf_gonad <- sapply(stringr::str_split(rownames(aggregated_ref), pattern="_", n=2), '[', 1) # Check that each cluster is present in all samples
alt_splitf_gonad <- sapply(stringr::str_split(rownames(aggregated_alt), pattern="_", n=2), '[', 1) # Check that each cluster is present in all samples

ref_pb_gonad <- split.data.frame(aggregated_ref, factor(ref_splitf_gonad)) %>% lapply(function(u) set_colnames(t(u), stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))
alt_pb_gonad <- split.data.frame(aggregated_alt, factor(alt_splitf_gonad)) %>% lapply(function(u) set_colnames(t(u), stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))


##############################################################################
########### Find MAR on all celltypes before splitting them up: ##############
##############################################################################

working_ac_gonad <- summed_ac_counts
working_ac_gonad <- working_ac_gonad %>%
  filter(ref_count >= 10 & alt_count >= 10)
working_ac_gonad$log10_ref <- log10(working_ac_gonad$ref_count)
working_ac_gonad$log10_alt <- log10(working_ac_gonad$alt_count)

working_ac_gonad <- working_ac_gonad %>%
  mutate(sample_id = case_when( 
    grepl("^M4", cells) ~ "M4",
    grepl("^M5", cells) ~ "M5",
    grepl("^M6", cells) ~ "M6",
    grepl("^F1", cells) ~ "F1",
    grepl("^F2", cells) ~ "F2",
    grepl("^F3", cells) ~ "F3"))

write.table(working_ac_gonad, "cellSNP_scAlleleCount/R_results/gonad_filtered10.txt", quote = F, row.names = F)
