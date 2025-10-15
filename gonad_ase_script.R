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

#write.table(working_ac_gonad, "scASE_results/gonad_filtered10.txt", quote = F, row.names = F)
write.table(working_ac_gonad, "cellSNP_scAlleleCount/R_results/gonad_filtered10.txt", quote = F, row.names = F)

plotA <- ggplot(working_ac_gonad, aes(x = sample_id, y = log10_ref)) + geom_violin() + ggtitle("Gonad") + ylab("log10 Ref Count")
plotB <- ggplot(working_ac_gonad, aes(x = sample_id, y = log10_alt)) + geom_violin() + ggtitle("Gonad") + ylab("log10 Alt Count")

plotA + plotB

#Distribution of the allele counts look equal across males vs. female samples

######################################
#### Filtering the allele counts #####
######################################

working_ac_gonad %>%
  subset(sample_id == "F1") %>% pull(ref_count) %>% quantile(probs = 0.2) # 16
working_ac_gonad %>%
  subset(sample_id == "F1") %>% pull(ref_count) %>% quantile(probs = 0.8) # 84
working_ac_gonad %>%
  subset(sample_id == "F1") %>% pull(alt_count) %>% quantile(probs = 0.2) # 12
working_ac_gonad %>%
  subset(sample_id == "F1") %>% pull(alt_count) %>% quantile(probs = 0.8) # 42


working_ac_gonad %>%
  subset(sample_id == "F2") %>% pull(ref_count) %>% quantile(probs = 0.2) # 16
working_ac_gonad %>%
  subset(sample_id == "F2") %>% pull(ref_count) %>% quantile(probs = 0.8) # 108
working_ac_gonad %>%
  subset(sample_id == "F2") %>% pull(alt_count) %>% quantile(probs = 0.2) # 13
working_ac_gonad %>%
  subset(sample_id == "F2") %>% pull(alt_count) %>% quantile(probs = 0.8) # 51

working_ac_gonad %>%
  subset(sample_id == "F3") %>% pull(ref_count) %>% quantile(probs = 0.2) # 16
working_ac_gonad %>%
  subset(sample_id == "F3") %>% pull(ref_count) %>% quantile(probs = 0.8) # 119
working_ac_gonad %>%
  subset(sample_id == "F3") %>% pull(alt_count) %>% quantile(probs = 0.2) # 13
working_ac_gonad %>%
  subset(sample_id == "F3") %>% pull(alt_count) %>% quantile(probs = 0.8) # 53

working_ac_gonad %>%
  subset(sample_id == "M4") %>% pull(ref_count) %>% quantile(probs = 0.2) # 14
working_ac_gonad %>%
  subset(sample_id == "M4") %>% pull(ref_count) %>% quantile(probs = 0.8) # 59
working_ac_gonad %>%
  subset(sample_id == "M4") %>% pull(alt_count) %>% quantile(probs = 0.2) # 12
working_ac_gonad %>%
  subset(sample_id == "M4") %>% pull(alt_count) %>% quantile(probs = 0.8) # 17

working_ac_gonad %>%
  subset(sample_id == "M5") %>% pull(ref_count) %>% quantile(probs = 0.2) # 13
working_ac_gonad %>%
  subset(sample_id == "M5") %>% pull(ref_count) %>% quantile(probs = 0.8) # 56
working_ac_gonad %>%
  subset(sample_id == "M5") %>% pull(alt_count) %>% quantile(probs = 0.2) # 12
working_ac_gonad %>%
  subset(sample_id == "M5") %>% pull(alt_count) %>% quantile(probs = 0.8) # 33

working_ac_gonad %>%
  subset(sample_id == "M6") %>% pull(ref_count) %>% quantile(probs = 0.2) # 14
working_ac_gonad %>%
  subset(sample_id == "M6") %>% pull(ref_count) %>% quantile(probs = 0.8) # 68
working_ac_gonad %>%
  subset(sample_id == "M6") %>% pull(alt_count) %>% quantile(probs = 0.2) # 12
working_ac_gonad %>%
  subset(sample_id == "M6") %>% pull(alt_count) %>% quantile(probs = 0.8) # 37

#working_ac_gonad_filtered <- working_ac_gonad

working_ac_gonad_filtered <- working_ac_gonad %>%
  filter(
      (sample_id == "F1" & ref_count < 75 & alt_count > 12) |
      (sample_id == "F2" & ref_count < 75 & alt_count > 12) |
      (sample_id == "F3" & ref_count < 75 & alt_count > 12) |
      (sample_id == "M4" & ref_count < 75 & alt_count > 12) |
      (sample_id == "M5" & ref_count < 75 & alt_count > 12) |
      (sample_id == "M6" & ref_count < 75 & alt_count > 12))

plotC <- ggplot(working_ac_gonad_filtered, aes(x = sample_id, y = log10_ref)) + geom_violin(trim = F) + ggtitle("Gonad Filtered") + ylab("log10 Ref Count")
plotD <- ggplot(working_ac_gonad_filtered, aes(x = sample_id, y = log10_alt)) + geom_violin(trim = F) + ggtitle("Gonad Filtered") + ylab("log10 Alt Count")

plotC + plotD


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

working_ac_gonad_filtered <- add_chrom_column(working_ac_gonad_filtered, gene_info)
working_ac_gonad_filtered$chrom_type <- ifelse(working_ac_gonad_filtered$chrom == "LG12", "Sex Chromosome", "Autosome")
working_ac_gonad_filtered$chrom_type[working_ac_gonad_filtered$geneid %in% gene_list] <- "Impacted Gene"

working_ac_gonad_A <- subset(working_ac_gonad_filtered, chrom_type == "Autosome")
working_ac_gonad_sc <- subset(working_ac_gonad_filtered, chrom_type == "Sex Chromosome")
working_ac_gonad_HI <- subset(working_ac_gonad_filtered, chrom_type == "Impacted Gene")

combo_working_ac_gonad <- rbind(working_ac_gonad_A, working_ac_gonad_HI, working_ac_gonad_sc)


celltype_ase <- combo_working_ac_gonad %>%
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

gonad_impactgenes <- subset(plot_ase, chrom_type == "Impacted Gene")
gonad_autosome <- subset(plot_ase, chrom_type == "Autosome")
gonad_sexchromo <- subset(plot_ase, chrom_type == "Sex Chromosome")

### Plotting:

gonad_impactgenes %>% filter(sex == "F") %>% pull(avg_mar) %>% median()
gonad_impactgenes %>% filter(sex == "M") %>% pull(avg_mar) %>% median()
gonad_autosome %>% filter(sex == "F") %>% pull(avg_mar) %>% median()
gonad_autosome %>% filter(sex == "M") %>% pull(avg_mar) %>% median()
gonad_sexchromo %>% filter(sex == "F") %>% pull(avg_mar) %>% median()
gonad_sexchromo %>% filter(sex == "M") %>% pull(avg_mar) %>% median()

plot_gonad_IG <- ggplot(gonad_impactgenes, aes(x = avg_mar)) + geom_density(aes(fill = sex), alpha = 0.5, adjust = 3) + 
  scale_fill_manual(values = c("firebrick3", "cornflowerblue")) +
  #scale_fill_manual(values = c("firebrick3", "hotpink", "pink", "navy", "blue", "cornflowerblue")) + 
  labs(fill = "Sex") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 12)) +
  ylab("Density") + xlab("Major Allele Ratio") + ggtitle("Gonad (Impacted Genes)") +
  geom_vline(aes(xintercept=0.5806452),color="firebrick4",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5493421),color="navy",linetype="dashed",size=0.8) +
  theme_classic() + theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    plot.title = element_text(size = 20, hjust = 0.5)
  )
plot_gonad_IG

plot_gonad_A <- ggplot(gonad_autosome, aes(x = avg_mar)) + geom_density(aes(fill = sex), alpha = 0.5) + 
  scale_fill_manual(values = c("firebrick3", "cornflowerblue")) +
  #scale_fill_manual(values = c("firebrick3", "hotpink", "pink", "navy", "blue", "cornflowerblue")) + 
  labs(fill = "Sex") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 12)) +
  ylab("Density") + xlab("Major Allele Ratio") + ggtitle("Gonad (Autosomes)") +
  geom_vline(aes(xintercept=0.5652174),color="firebrick4",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5636364),color="navy",linetype="dashed",size=0.8) +
  theme_classic() + theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    plot.title = element_text(size = 20, hjust = 0.5)
  )
plot_gonad_A

plot_gonad_SC <- ggplot(gonad_sexchromo, aes(x = avg_mar)) + geom_density(aes(fill = sex), alpha = 0.5) + 
  scale_fill_manual(values = c("firebrick3", "cornflowerblue")) +
  #scale_fill_manual(values = c("firebrick3", "hotpink", "pink", "navy", "blue", "cornflowerblue")) + 
  labs(fill = "Sex") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 12)) +
  ylab("Density") + xlab("Major Allele Ratio") + ggtitle("Gonad (Sex Chromosome)") +
  geom_vline(aes(xintercept=0.5681177),color="firebrick4",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5680604),color="navy",linetype="dashed",size=0.8) +
  theme_classic() + theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    plot.title = element_text(size = 20, hjust = 0.5)
  )
plot_gonad_SC

gonad_plots <- (plot_gonad_IG | plot_gonad_A | plot_gonad_SC) + plot_layout(guides = "collect")
gonad_plots

wilcox.test(gonad_impactgenes$avg_mar ~ gonad_impactgenes$sex, p.adjust.methods = "bonferroni")
wilcox.test(gonad_autosome$avg_mar ~ gonad_autosome$sex, p.adjust.methods = "bonferroni")
wilcox.test(gonad_sexchromo$avg_mar ~ gonad_sexchromo$sex, p.adjust.methods = "bonferroni")



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
write.csv(res_df_male, "scASE_results/gonad_male_all.csv", row.names = T)
write.csv(sig_res_df_male, "scASE_results/gonad_male_all_sig.csv", row.names = T)

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
write.csv(res_df_female, "scASE_results/gonad_female_all.csv", row.names = T)
write.csv(sig_res_df_female, "scASE_results/gonad_female_all_sig.csv", row.names = T)


########### To run this for all celltypes in your matrix ##############

# Write a function so you can pull out every matrix from your ^_pb list:
# Separate out male vs. female

run_scASE_by_sex <- function(cell_type, ref_list, alt_list, output_dir = "scASE_results/gonad/") {
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

cell_types <- names(ref_pb_gonad)

results_by_sex <- lapply(cell_types, function(ct) {
  run_scASE_by_sex(ct, ref_pb_gonad, alt_pb_gonad)
})

# Since "Germ" does not have the respective "M", this function stops running after Germ_F
# Redo the "cell_types" to continue running the other cell_types

cell_types_no_germ <- c("Granulosa", "Hemoglobin rich", "Macrophage", "Monocyte", "Spermatid", "Spermatocytes")

results_by_sex_no_germ <- lapply(cell_types_no_germ, function(ct) {
  run_scASE_by_sex(ct, ref_pb_gonad, alt_pb_gonad)
})

cell_types_no_macro <- c("Monocyte", "Spermatid", "Spermatocytes")

results_by_sex_no_macro <- lapply(cell_types_no_macro, function(ct) {
  run_scASE_by_sex(ct, ref_pb_gonad, alt_pb_gonad)
})


#names(results_by_sex) <- cell_types

##########################################
#### Macrophage did not run properly: ####
##########################################

ref_mat <- as.matrix(ref_pb_gonad[["Macrophage"]])
alt_mat <- as.matrix(alt_pb_gonad[["Macrophage"]])

ref_df <- as.data.frame(ref_mat) %>%
  tibble::rownames_to_column("gene_id") %>%
  pivot_longer(-gene_id, names_to = "cell", values_to = "ref")

alt_df <- as.data.frame(alt_mat) %>%
  tibble::rownames_to_column("gene_id") %>%
  pivot_longer(-gene_id, names_to = "cell", values_to = "alt")

combined_df <- ref_df %>%
  mutate(alt = alt_df$alt,
         sex = ifelse(grepl("F", cell), "F", "M")) 

long_df <- combined_df %>%
  pivot_longer(cols = c(ref, alt), names_to = "allele", values_to = "count") %>%
  mutate(
    cell = factor(cell),
    allele = factor(allele, levels = c("ref", "alt")),
    sex = factor(sex)
  ) %>%
  unite(sample, cell, allele, sep = "_", remove = FALSE) %>%
  na.omit()

sex_F <- filter(long_df, sex == "F")
count_mat <- sex_F %>%
  select(gene_id, sample, count) %>%
  pivot_wider(names_from = sample, values_from = count, values_fill = 0) %>%
  column_to_rownames("gene_id") %>%
  as.matrix()
coldata <- sex_F %>%
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

write.csv(res_df, "scASE_results/gonad/Macrophage_F.csv", row.names = T)
write.csv(sig_res_df, "scASE_results/gonad/Macrophage_sig_F.csv", row.names = T)

  # Run for Male cells:

sex_M <- filter(long_df, sex == "M")
count_mat <- sex_M %>%
  select(gene_id, sample, count) %>%
  pivot_wider(names_from = sample, values_from = count, values_fill = 0) %>%
  column_to_rownames("gene_id") %>%
  as.matrix()
#Given the number of 0s per gene
#Error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  : 
# every gene contains at least one zero, cannot compute log geometric means
#Add a "1" to all the data cells
count_mat <- count_mat+1
coldata <- sex_M %>%
  distinct(sample, allele, cell) %>%
  as.data.frame()
rownames(coldata) <- coldata$sample

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = coldata,
  design = ~ allele
)

dds <- DESeq(estimateSizeFactors(dds, type = 'iterate')) #To account for the added "1"
res <- results(dds, contrast = c("allele", "alt", "ref"))
res_df <- as.data.frame(res) %>%
  filter(log2FoldChange != 0)
sig_res_df <- as.data.frame(res) %>%
  filter(log2FoldChange != 0, padj < 0.05)

write.csv(res_df, "scASE_results/gonad/Macrophage_M.csv", row.names = T)
write.csv(sig_res_df, "scASE_results/gonad/Macrophage_sig_M.csv", row.names = T)



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


    # 1. Endocrine

endo_F <- read.csv("scASE_results/gonad/Endocrine_F.csv")
endo_M <- read.csv("scASE_results/gonad/Endocrine_M.csv")
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
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.5643516
endo_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.5472773
endo_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.5708093


endo_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.5919882
endo_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.576966
endo_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.5975942

write.csv(endo_combo_M, "scASE_results/gonad/split_chromo/Endocrine_M.csv", quote = FALSE, row.names = FALSE)
write.csv(endo_combo_F, "scASE_results/gonad/split_chromo/Endocrine_F.csv", quote = FALSE, row.names = FALSE)


plot_endo_M <- ggplot(endo_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 2) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 8)) +
  geom_vline(aes(xintercept=0.5643516),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5472773),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6047359),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Gonad (Endocrine M)")
plot_endo_M

plot_endo_F <- ggplot(endo_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.5) + 
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 8)) +
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  geom_vline(aes(xintercept=0.5919882),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.576966),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5975942),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Gonad (Endocrine F)")
plot_endo_F


    # 2. Germ

germ_F <- read.csv("scASE_results/gonad/Germ_F.csv")
germ_F$sex = "F"


germ_F$prop <- (2^germ_F$log2FoldChange)/(1+2^germ_F$log2FoldChange)
germ_F$mar <- pmax(germ_F$prop, 1 - germ_F$prop)
germ_F <- add_chrom_column(germ_F, gene_info)
germ_F$chrom_type <- ifelse(germ_F$chrom == "LG12", "Sex Chromosome", "Autosome")
germ_F$chrom_type[germ_F$X %in% gene_list] <- "Impacted Gene"

germ_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.6017565
germ_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.5631761
germ_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) #  0.6219763

write.csv(germ_F, "scASE_results/gonad/split_chromo/Germ_F.csv", quote = FALSE, row.names = FALSE)

plot_germ_F <- ggplot(germ_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.5) + 
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 8)) +
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  geom_vline(aes(xintercept=0.6017565),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5631761),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6219763),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Gonad (Germ cells F)")
plot_germ_F


    # 3. Granulosa

granu_F <- read.csv("scASE_results/gonad/Granulosa_F.csv")
granu_M <- read.csv("scASE_results/gonad/Granulosa_M.csv")
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
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.6053268
granu_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.7324507
granu_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6525655


granu_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.5967189
granu_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.5777726
granu_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.5871192

write.csv(granu_combo_M, "scASE_results/gonad/split_chromo/Granulosa_M.csv", quote = FALSE, row.names = FALSE)
write.csv(granu_combo_F, "scASE_results/gonad/split_chromo/Granulosa_F.csv", quote = FALSE, row.names = FALSE)


plot_granu_M <- ggplot(granu_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.5) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 8)) +
  geom_vline(aes(xintercept=0.6053268),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.7324507),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6525655),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Gonad (Granulosa M)")
plot_granu_M

plot_granu_F <- ggplot(granu_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1) + 
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 8)) +
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  geom_vline(aes(xintercept=0.5967189),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5777726),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5871192),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Gonad (Granulosa F)")
plot_granu_F


    # 4. Hemoglobin Rich

hemo_F <- read.csv("scASE_results/gonad/Hemoglobin rich_F.csv")
hemo_M <- read.csv("scASE_results/gonad/Hemoglobin rich_M.csv")
hemo_F$sex = "F"
hemo_M$sex = "M"
hemo_combo <- rbind(hemo_F, hemo_M)

hemo_combo$prop <- (2^hemo_combo$log2FoldChange)/(1+2^hemo_combo$log2FoldChange)
hemo_combo$mar <- pmax(hemo_combo$prop, 1 - hemo_combo$prop)
hemo_combo <- add_chrom_column(hemo_combo, gene_info)
hemo_combo$chrom_type <- ifelse(hemo_combo$chrom == "LG12", "Sex Chromosome", "Autosome")
hemo_combo$chrom_type[hemo_combo$X %in% gene_list] <- "Impacted Gene"
hemo_combo_M <- subset(hemo_combo, sex == "M")
hemo_combo_F <- subset(hemo_combo, sex == "F")

hemo_combo_M %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.610741
hemo_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.5656188
hemo_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6060262


hemo_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.6284994
hemo_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.6384447
hemo_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6550843

write.csv(hemo_combo_M, "scASE_results/gonad/split_chromo/Hemoglobin_M.csv", quote = FALSE, row.names = FALSE)
write.csv(hemo_combo_F, "scASE_results/gonad/split_chromo/Hemoglobin_F.csv", quote = FALSE, row.names = FALSE)


plot_hemo_M <- ggplot(hemo_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 8)) +
  geom_vline(aes(xintercept=0.610741),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5656188),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6060262),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Gonad (Hemoglobin rich cells M)")
plot_hemo_M

plot_hemo_F <- ggplot(hemo_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.2) + 
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 8)) +
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  geom_vline(aes(xintercept=0.6284994),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6384447),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6550843),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Gonad (Hemoglobin rich cells F)")
plot_hemo_F


    # 5. Macrophage
# Not enough gene counts for the male macrophage cells to be meaningful

macro_F <- read.csv("scASE_results/gonad/Macrophage_F.csv")
macro_M <- read.csv("scASE_results/gonad/Macrophage_M.csv")
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
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.5105408
macro_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.5105408
macro_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.5105408

macro_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.587665
macro_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.5751406
macro_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6111177

write.csv(macro_combo_M, "scASE_results/gonad/split_chromo/Macrophage_M.csv", quote = FALSE, row.names = FALSE)
write.csv(macro_combo_F, "scASE_results/gonad/split_chromo/Macrophage_F.csv", quote = FALSE, row.names = FALSE)


plot_macro_M <- ggplot(macro_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 4.5) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 8)) +
  geom_vline(aes(xintercept=0.5105408),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5105408),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5105408),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Gonad (Macrophage M)")
plot_macro_M

plot_macro_F <- ggplot(macro_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.5) + 
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 8)) +
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  geom_vline(aes(xintercept=0.587665),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5751406),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6111177),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Gonad (Macrophage F)")
plot_macro_F


    # 6. Monocyte
mono_F <- read.csv("scASE_results/gonad/Monocyte_F.csv")
mono_M <- read.csv("scASE_results/gonad/Monocyte_M.csv")
mono_F$sex = "F"
mono_M$sex = "M"
mono_combo <- rbind(mono_F, mono_M)

mono_combo$prop <- (2^mono_combo$log2FoldChange)/(1+2^mono_combo$log2FoldChange)
mono_combo$mar <- pmax(mono_combo$prop, 1 - mono_combo$prop)
mono_combo <- add_chrom_column(mono_combo, gene_info)
mono_combo$chrom_type <- ifelse(mono_combo$chrom == "LG12", "Sex Chromosome", "Autosome")
mono_combo$chrom_type[mono_combo$X %in% gene_list] <- "Impacted Gene"
mono_combo_M <- subset(mono_combo, sex == "M")
mono_combo_F <- subset(mono_combo, sex == "F")

mono_combo_M %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.5911961
mono_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) #  0.5003621
mono_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.5695795

mono_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.6160379
mono_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.5800464
mono_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6217076

write.csv(mono_combo_M, "scASE_results/gonad/split_chromo/Monocyte_M.csv", quote = FALSE, row.names = FALSE)
write.csv(mono_combo_F, "scASE_results/gonad/split_chromo/Monocyte_F.csv", quote = FALSE, row.names = FALSE)


plot_mono_M <- ggplot(mono_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.5) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 8)) +
  geom_vline(aes(xintercept=0.5911961),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5003621),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5695795),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Gonad (Monocyte M)")
plot_mono_M

plot_mono_F <- ggplot(mono_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.5) + 
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 8)) +
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  geom_vline(aes(xintercept=0.6160379),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5800464),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6217076),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Gonad (Monocyte F)")
plot_mono_F


    # 7. spermatid

spermatid_F <- read.csv("scASE_results/gonad/Spermatid_F.csv")
spermatid_M <- read.csv("scASE_results/gonad/Spermatid_M.csv")
spermatid_F$sex = "F"
spermatid_M$sex = "M"
spermatid_combo <- rbind(spermatid_F, spermatid_M)

spermatid_combo$prop <- (2^spermatid_combo$log2FoldChange)/(1+2^spermatid_combo$log2FoldChange)
spermatid_combo$mar <- pmax(spermatid_combo$prop, 1 - spermatid_combo$prop)
spermatid_combo <- add_chrom_column(spermatid_combo, gene_info)
spermatid_combo$chrom_type <- ifelse(spermatid_combo$chrom == "LG12", "Sex Chromosome", "Autosome")
spermatid_combo$chrom_type[spermatid_combo$X %in% gene_list] <- "Impacted Gene"
spermatid_combo_M <- subset(spermatid_combo, sex == "M")
spermatid_combo_F <- subset(spermatid_combo, sex == "F")

spermatid_combo_M %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.5517695
spermatid_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) #  0.5688501
spermatid_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) #  0.5595575

spermatid_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.6339354
spermatid_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.6215608
spermatid_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6339354

write.csv(spermatid_combo_M, "scASE_results/gonad/split_chromo/Spermatid_M.csv", quote = FALSE, row.names = FALSE)
write.csv(spermatid_combo_F, "scASE_results/gonad/split_chromo/Spermatid_F.csv", quote = FALSE, row.names = FALSE)


plot_spermatid_M <- ggplot(spermatid_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 3.5) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 8)) +
  geom_vline(aes(xintercept=0.5517695),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5688501),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5595575),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Gonad (Spermatid M)")
plot_spermatid_M

plot_spermatid_F <- ggplot(spermatid_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.5) + 
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 8)) +
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  geom_vline(aes(xintercept=0.6339354),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6215608),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6339354),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Gonad (Spermatid F)")
plot_spermatid_F

    # 8. spermatocytes


spermatocyte_F <- read.csv("scASE_results/gonad/Spermatocytes_F.csv")
spermatocyte_M <- read.csv("scASE_results/gonad/Spermatocytes_M.csv")
spermatocyte_F$sex = "F"
spermatocyte_M$sex = "M"
spermatocyte_combo <- rbind(spermatocyte_F, spermatocyte_M)

spermatocyte_combo$prop <- (2^spermatocyte_combo$log2FoldChange)/(1+2^spermatocyte_combo$log2FoldChange)
spermatocyte_combo$mar <- pmax(spermatocyte_combo$prop, 1 - spermatocyte_combo$prop)
spermatocyte_combo <- add_chrom_column(spermatocyte_combo, gene_info)
spermatocyte_combo$chrom_type <- ifelse(spermatocyte_combo$chrom == "LG12", "Sex Chromosome", "Autosome")
spermatocyte_combo$chrom_type[spermatocyte_combo$X %in% gene_list] <- "Impacted Gene"
spermatocyte_combo_M <- subset(spermatocyte_combo, sex == "M")
spermatocyte_combo_F <- subset(spermatocyte_combo, sex == "F")

spermatocyte_combo_M %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) #  0.5649256
spermatocyte_combo_M %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) #  0.5478292
spermatocyte_combo_M %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) #  0.5709225

spermatocyte_combo_F %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar)) # 0.6210637
spermatocyte_combo_F %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar)) # 0.8686594
spermatocyte_combo_F %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar)) # 0.6771388

write.csv(spermatocyte_combo_M, "scASE_results/gonad/split_chromo/Spermatocytes_M.csv", quote = FALSE, row.names = FALSE)
write.csv(spermatocyte_combo_F, "scASE_results/gonad/split_chromo/Spermatocytes_F.csv", quote = FALSE, row.names = FALSE)


plot_spermatocyte_M <- ggplot(spermatocyte_combo_M, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 3) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 8)) +
  geom_vline(aes(xintercept=0.5649256),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5478292),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.5709225),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Gonad (Spermatocytes M)")
plot_spermatocyte_M

plot_spermatocyte_F <- ggplot(spermatocyte_combo_F, aes(x = mar)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.2) + 
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 8)) +
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  geom_vline(aes(xintercept=0.6210637),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.8686594),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.6771388),color="darkorchid4",linetype="dashed",size=0.8) + 
  theme_classic() + ggtitle("Gonad (Spermatocytes F)")
plot_spermatocyte_F


######################################################

#### Combo your plots ####

library(patchwork)

gonad_celltype_plots_M <- wrap_plots(list(plot_endo_M, plot_hemo_M, plot_mono_M, plot_spermatid_M, plot_spermatocyte_M), 
                                     ncol = 2, nrow = 3) + plot_layout(guides = "collect")

gonad_celltype_plots_M

gonad_celltype_plots_F <- wrap_plots(list(plot_endo_F, plot_hemo_F, plot_mono_F, plot_granu_F, plot_macro_F, plot_germ_F), 
                                     ncol = 2, nrow = 3) +  plot_layout(guides = "collect")

gonad_celltype_plots_F


######################################################
######## Combine celltypes to test for ASE ###########
######################################################


all_celltype_M <- rbind(endo_combo_M, hemo_combo_M, mono_combo_M,
                        spermatid_combo_M, spermatocyte_combo_M)

all_celltype_F <- rbind(endo_combo_F, hemo_combo_F, mono_combo_F,
                        macro_combo_F, granu_combo_F, germ_F)


all_celltype_M_means <- all_celltype_M %>%
  group_by(X, chrom_type) %>%
  summarise(mar_mean = mean(mar), .groups = "drop")

all_celltype_F_means <- all_celltype_F %>%
  group_by(X, chrom_type) %>%
  summarise(mar_mean = mean(mar), .groups = "drop")



all_celltype_M_means %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar_mean)) #  0.569
all_celltype_M_means %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar_mean)) # 0.564
all_celltype_M_means %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar_mean)) # 0.574

all_celltype_F_means %>% 
  filter(chrom_type == "Autosome") %>% summarise(median(mar_mean)) # 0.613
all_celltype_F_means %>% 
  filter(chrom_type == "Impacted Gene") %>% summarise(median(mar_mean)) # 0.594
all_celltype_F_means %>% 
  filter(chrom_type == "Sex Chromosome") %>% summarise(median(mar_mean)) # 0.613


plot_celltype_M <- ggplot(all_celltype_M_means, aes(x = mar_mean)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.5) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 15)) +
  geom_vline(aes(xintercept=0.569),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.564),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept= 0.574),color="darkorchid4",linetype="dashed",size=0.8) +
  theme_classic() + ggtitle("Gonad M (all celltypes)")
plot_celltype_M


plot_celltype_F <- ggplot(all_celltype_F_means, aes(x = mar_mean)) + geom_density(aes(fill = chrom_type), alpha = 0.5, adjust = 1.5) + 
  scale_fill_manual(values = c("grey40", "orange", "darkorchid4")) + labs(fill = "Chromosome Type") +
  scale_x_continuous(lim = c(0.49, 1.0), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 15)) +
  geom_vline(aes(xintercept=0.613),color="grey20",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept=0.594),color="orange",linetype="dashed",size=0.8) +
  geom_vline(aes(xintercept= 0.613),color="darkorchid4",linetype="dashed",size=0.8) +
  theme_classic() + ggtitle("Gonad F (all celltypes)")
plot_celltype_F


combo_everything <- wrap_plots(plot_celltype_M, plot_celltype_F) +
  plot_layout(guides = "collect")

combo_everything

