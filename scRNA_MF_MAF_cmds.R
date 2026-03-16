######################################################
############ Prepping Iulia's scRNA seq ##############
############### For M:F MAF Plots ####################
######################################################

setwd("//files.zoology.ubc.ca/ljmfong/flex")
getwd()

rm(list=ls())
ls() 


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

######################################################
#####################    HEART    ####################
######################################################

load("iulia_scRNA/heart_umap_3_rename_clusters_dbrmv.RData")

counts_heart <- heart_umap_3_rename_clusters_dbrmv@assays$RNA@counts
metadata_heart <- heart_umap_3_rename_clusters_dbrmv@meta.data
metadata_heart$CellType <- factor(heart_umap_3_rename_clusters_dbrmv@active.ident)

cleaned_df <- read.table(file = "cellSNP_scAlleleCount/R_results/snpfilt_heart_cleaned_df_gene_attached.txt", header = T)

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
  group_by(cells, geneid, CHROM, POS) %>%
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
cleaned_group_name_ac_heart$summed_counts <- (cleaned_group_name_ac_heart$ref_count) + (cleaned_group_name_ac_heart$alt_count)

working_ac_heart <- cleaned_group_name_ac_heart %>%
  filter(summed_counts >= 10)

working_ac_heart <- working_ac_heart %>%
  mutate(sample_id = case_when( 
    grepl("^M1", cells) ~ "M1",
    grepl("^M2", cells) ~ "M2",
    grepl("^M3", cells) ~ "M3",
    grepl("^F1", cells) ~ "F1",
    grepl("^F2", cells) ~ "F2",
    grepl("^F3", cells) ~ "F3"))

working_ac_heart$tissue <- "heart"


######################################################
#####################    LIVER    ####################
######################################################

load("iulia_scRNA/liver_umap_3_rename_clusters_dbrmv.RData")

#### 3. Get differential gene expression for separate cell types ####
counts_liver <- liver_umap_3_rename_clusters_dbrmv@assays$RNA@counts
metadata_liver <- liver_umap_3_rename_clusters_dbrmv@meta.data
metadata_liver$CellType <- factor(liver_umap_3_rename_clusters_dbrmv@active.ident)

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
  group_by(cells, geneid, CHROM, POS) %>%
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
cleaned_group_name_ac_liver$summed_counts <- (cleaned_group_name_ac_liver$ref_count) + (cleaned_group_name_ac_liver$alt_count)

working_ac_liver <- cleaned_group_name_ac_liver %>%
  filter(summed_counts >= 10)

working_ac_liver <- working_ac_liver %>%
  mutate(sample_id = case_when( 
    grepl("^M1", cells) ~ "M1",
    grepl("^M2", cells) ~ "M2",
    grepl("^M3", cells) ~ "M3",
    grepl("^F1", cells) ~ "F1",
    grepl("^F2", cells) ~ "F2",
    grepl("^F3", cells) ~ "F3"))

working_ac_liver$tissue <- "liver"


######################################################
#####################    SKIN    #####################
######################################################

load("iulia_scRNA/skin_umap_3_rename_clusters_dbrmv.RData")

counts_skin <- skin_umap_3_rename_clusters_dbrmv@assays$RNA@counts
metadata_skin <- skin_umap_3_rename_clusters_dbrmv@meta.data
metadata_skin$CellType <- factor(skin_umap_3_rename_clusters_dbrmv@active.ident)

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


summed_ac_counts <- cleaned_df %>%
  group_by(cells, geneid, CHROM, POS) %>%
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
cleaned_group_name_ac_skin$summed_counts <- (cleaned_group_name_ac_skin$ref_count) + (cleaned_group_name_ac_skin$alt_count)

working_ac_skin <- cleaned_group_name_ac_skin %>%
  filter(summed_counts >= 10)

working_ac_skin <- working_ac_skin %>%
  mutate(sample_id = case_when( 
    grepl("^M1", cells) ~ "M1",
    grepl("^M2", cells) ~ "M2",
    grepl("^M3", cells) ~ "M3",
    grepl("^F1", cells) ~ "F1",
    grepl("^F2", cells) ~ "F2",
    grepl("^F3", cells) ~ "F3"))

working_ac_skin$tissue <- "skin"


######################################################
####################    GONAD    #####################
######################################################

load("iulia_scRNA/gonad_umap_3_rename_clusters.RData")

counts_gonad <- gonad_umap_3_rename_clusters@assays$RNA@counts
metadata_gonad <- gonad_umap_3_rename_clusters@meta.data
metadata_gonad$CellType <- factor(gonad_umap_3_rename_clusters@active.ident)

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

summed_ac_counts <- cleaned_df %>%
  group_by(cells, CHROM, POS, geneid) %>%
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
cleaned_group_name_ac_gonad$summed_counts <- (cleaned_group_name_ac_gonad$ref_count) + (cleaned_group_name_ac_gonad$alt_count)

working_ac_gonad <- cleaned_group_name_ac_gonad %>%
  filter(summed_counts >= 10)

working_ac_gonad <- working_ac_gonad %>%
  mutate(sample_id = case_when( 
    grepl("^M4", cells) ~ "M4",
    grepl("^M5", cells) ~ "M5",
    grepl("^M6", cells) ~ "M6",
    grepl("^F1", cells) ~ "F1",
    grepl("^F2", cells) ~ "F2",
    grepl("^F3", cells) ~ "F3"))

working_ac_gonad$tissue <- "gonad"


#################################################
############### M:F MAF Windows #################
#################################################

library(zoo)
library(RcppRoll)

movingaverage <- function (x, window) {
  ma <- roll_mean(x, window, fill = NA)
}
windowsize <- 10 

  #### Prepping file #####

chrm_size <- read.table(file.choose(), sep = "\t") #headers are = V1, V2, V3
# File is: /home/ljmfong/guppy_chrms_len_10kb.txt

somatic_all <- rbind(working_ac_heart, working_ac_liver, working_ac_skin)
gonad_all <- working_ac_gonad
all_tissue <- rbind(somatic_all, working_ac_gonad)

somatic_allele_counts <- somatic_all %>%  filter(ref_count <= 49 & alt_count <= 104)
gonad_allele_counts <- gonad_all %>%  filter(ref_count <= 41 & alt_count <= 78)

combo_all_filt <- rbind(somatic_allele_counts, gonad_allele_counts)
combo_all_filt$sex <- ifelse(grepl("^F", combo_all_filt$sample_id), "F", "M")
somatic_tiss <- subset(combo_all_filt, tissue != "gonad")
somatic_tiss <- somatic_tiss %>% group_by(geneid) %>%
  filter(sum(sex == "F") >= 3, sum(sex == "M") >= 3) %>% ungroup()

gonad_tiss <- subset(combo_all_filt, tissue == "gonad")
gonad_tiss <- gonad_tiss %>% group_by(geneid) %>%
  filter(sum(sex == "F") >= 3, sum(sex == "M") >= 3) %>% ungroup()


combo_all_filt <- rbind(somatic_tiss, gonad_tiss)
combo_all_filt$total_counts <- combo_all_filt$ref_count + combo_all_filt$alt_count
combo_all_filt$mar <- pmax(combo_all_filt$ref_count, combo_all_filt$alt_count)/(combo_all_filt$total_counts)
combo_all_filt <- na.omit(combo_all_filt)

somatic_tiss_w_mar <- subset(combo_all_filt, tissue != "gonad")
gonad_tiss_w_mar <- subset(combo_all_filt, tissue == "gonad")

fem_somatic_mar <- subset(somatic_tiss_w_mar, sex == "F")
fem_gonad_mar <- subset(gonad_tiss_w_mar, sex == "F")
male_somatic_mar <- subset(somatic_tiss_w_mar, sex == "M")
male_gonad_mar <- subset(gonad_tiss_w_mar, sex == "M")

  ##### Put the tissues into windows #####

fem_som_sorted <- fem_somatic_mar[order(fem_somatic_mar$CHROM, fem_somatic_mar$POS),]
fem_gond_sorted <- fem_gonad_mar[order(fem_gonad_mar$CHROM, fem_gonad_mar$POS),]
male_som_sorted <- male_somatic_mar[order(male_somatic_mar$CHROM, male_somatic_mar$POS),]
male_gond_sorted <- male_gonad_mar[order(male_gonad_mar$CHROM, male_gonad_mar$POS),]

fem_som_smooth <- fem_som_sorted %>%
  mutate(window = floor(POS/100000)*100000) %>% group_by(CHROM, window) %>%
  summarise(mean_mar_100kb = mean(mar, na.rm = T), .groups = "drop")
male_som_smooth <- male_som_sorted %>%
  mutate(window = floor(POS/100000)*100000) %>% group_by(CHROM, window) %>%
  summarise(mean_mar_100kb = mean(mar, na.rm = T), .groups = "drop")
 
matching_genes_som <- left_join(fem_som_smooth, male_som_smooth, by=c("CHROM","window"))

fem_gond_smooth <- fem_gond_sorted %>%
  mutate(window = floor(POS/100000)*100000) %>% group_by(CHROM, window) %>%
  summarise(mean_mar_100kb = mean(mar, na.rm = T), .groups = "drop")
male_gond_smooth <- male_gond_sorted %>%
  mutate(window = floor(POS/100000)*100000) %>% group_by(CHROM, window) %>%
  summarise(mean_mar_100kb = mean(mar, na.rm = T), .groups = "drop")

matching_genes_gond <- left_join(fem_gond_smooth, male_gond_smooth, by=c("CHROM","window"))


  #### Running MF MAF ####

matching_genes_som$MF_MAF <- log2(matching_genes_som$mean_mar_100kb.y/matching_genes_som$mean_mar_100kb.x)
matching_genes_gond$MF_MAF <- log2(matching_genes_gond$mean_mar_100kb.y/matching_genes_gond$mean_mar_100kb.x)

matching_genes_som <- na.omit(matching_genes_som)
matching_genes_gond <- na.omit(matching_genes_gond)

somatic_mar_auto <- subset(matching_genes_som, CHROM != "LG12")
gonad_mar_auto <- subset(matching_genes_gond, CHROM != "LG12")
somatic_mar_sexchromo <- subset(matching_genes_som, CHROM == "LG12")
gonad_mar_sexchromo <- subset(matching_genes_gond, CHROM == "LG12")


  #### Permute for CI ####

MFautopermute_som <- replicate(1000,mean(sample(matching_genes_som$MF_MAF,windowsize,replace = F)))
MFautoI25cov_som <- quantile(MFautopermute_som, c(.025, .5, .975))[[1]]
MFautoI25cov_som

MFautoI975cov_som <- quantile(MFautopermute_som, c(.025, .5, .975))[[3]]
MFautoI975cov_som

MFautopermute_gond <- replicate(1000,mean(sample(matching_genes_gond$MF_MAF,windowsize,replace = F)))
MFautoI25cov_gond <- quantile(MFautopermute_gond, c(.025, .5, .975))[[1]]
MFautoI25cov_gond

MFautoI975cov_gond <- quantile(MFautopermute_gond, c(.025, .5, .975))[[3]]
MFautoI975cov_gond


  ##### Soma ######

smoothline_som <- movingaverage(somatic_mar_sexchromo$MF_MAF, windowsize)
sort_som_chrm12 <- somatic_mar_sexchromo[order(somatic_mar_sexchromo$window),]
line_df <- as.data.frame(smoothline_som)
line_df$start <- sort_som_chrm12$window
line_df <- na.omit(line_df)

ggplot(sort_som_chrm12, aes(x= window/1000000, y= MF_MAF)) +
  geom_rect(xmax=27,xmin=-10,ymin=MFautoI25cov_som,ymax=MFautoI975cov_som, fill="grey70", alpha=0.08)+
  geom_point(colour="grey50",fill="grey50", alpha=0.5, cex=1.5) +
  geom_line(data = line_df, aes(x=start/1000000, y=smoothline_som), colour = "black", size=1.5) +
  coord_cartesian(ylim=c(-2,2)) +
  coord_cartesian(xlim=c(0,27)) +
  theme(panel.grid.minor = element_blank(), panel.background = element_blank())+
  theme(
    text=element_text(size=12),
    plot.margin = unit(c(1.5,1.5,1.5,1.5),"lines"),
    axis.line.y = element_line(color="black", size = 0.3),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.text.x = element_text(size=20),
    axis.text.y = element_text(size=20), axis.title = element_text(size = 20)
  ) +
  ylim(-1.2,1.2) + 
  scale_x_continuous(expand = c(0,0), breaks=seq(0,27,5)) +
  xlab('Chromosome 12 (Mb)') +
  ylab(expression('M:F log'[2]*' Major Allele Ratio'))


  ##### Gonad #####

smoothline_gond <- movingaverage(gonad_mar_sexchromo$MF_MAF, windowsize)
sort_gond_chrm12 <- gonad_mar_sexchromo[order(gonad_mar_sexchromo$window),]
line_df_g <- as.data.frame(smoothline_gond)
line_df_g$start <- sort_gond_chrm12$window
line_df_g <- na.omit(line_df_g)

ggplot(sort_gond_chrm12, aes(x= window/1000000, y= MF_MAF)) +
  geom_rect(xmax=27,xmin=-10,ymin=MFautoI25cov_gond,ymax=MFautoI975cov_gond, fill="grey70", alpha=0.08)+
  geom_point(colour="grey50",fill="grey50", alpha=0.5, cex=1.5) +
  geom_line(data = line_df_g, aes(x=start/1000000, y=smoothline_gond), colour = "black", size=1.5) +
  coord_cartesian(ylim=c(-2,2)) +
  coord_cartesian(xlim=c(0,27)) +
  theme(panel.grid.minor = element_blank(), panel.background = element_blank())+
  theme(
    text=element_text(size=12),
    plot.margin = unit(c(1.5,1.5,1.5,1.5),"lines"),
    axis.line.y = element_line(color="black", size = 0.3),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.text.x = element_text(size=20),
    axis.text.y = element_text(size=20), axis.title = element_text(size = 20)
  ) +
  ylim(-1.2,1.2) + 
  scale_x_continuous(expand = c(0,0), breaks=seq(0,27,5)) +
  xlab('Chromosome 12 (Mb)') +
  ylab(expression('M:F log'[2]*' Major Allele Ratio'))

