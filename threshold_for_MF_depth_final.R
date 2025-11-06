#### Plotting the M:F Read Depth for genes in P. reticulata ####

rm(list=ls())
ls() 

library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)

setwd("\\\\files.zoology.ubc.ca/ljmfong/flex/new_Ret_SNPs")

####### 1. Read in the file: #######
####### Look at snps first #########

#Remove genes that are putative duplicates a.k.a those that you are interested in through HI snps

gene_list_autosome <- c("ENSPREG00000000149", "ENSPREG00000000425", "ENSPREG00000000486", "ENSPREG00000020299", 
               "ENSPREG00000021037", "ENSPREG00000014816", "ENSPREG00000005754", "ENSPREG00000013757", "ENSPREG00000003677",
               "ENSPREG00000015247", "ENSPREG00000006517", "ENSPREG00000011757", "ENSPREG00000013879", "ENSPREG00000004959",
               "ENSPREG00000021017", "ENSPREG00000001126", "ENSPREG00000007572", "ENSPREG00000003854", "ENSPREG00000019543",
               "ENSPREG00000003921", "ENSPREG00000012441", "ENSPREG00000012471", "ENSPREG00000004135", "ENSPREG00000017067")
gene_list_sexchromo <- c("ENSPREG00000019446", "ENSPREG00000006420", "ENSPREG00000018782", "ENSPREG00000020252", "ENSPREG00000019378",
                         "ENSPREG00000019435", "ENSPREG00000019478", "ENSPREG00000019540", "ENSPREG00000019673",
                         "ENSPREG00000019774", "ENSPREG00000019794", "ENSPREG00000019838", "ENSPREG00000019848", "ENSPREG00000019881",
                         "ENSPREG00000019905", "ENSPREG00000019944", "ENSPREG00000019970", "ENSPREG00000019988",
                         "ENSPREG00000020032", "ENSPREG00000020049", "ENSPREG00000020126", "ENSPREG00000020131", "ENSPREG00000020223",
                         "ENSPREG00000020570", "ENSPREG00000001200", "ENSPREG00000001330", "ENSPREG00000001446", "ENSPREG00000001498", "ENSPREG00000001901", "ENSPREG00000005501",
                         "ENSPREG00000005616", "ENSPREG00000005759", "ENSPREG00000005809", "ENSPREG00000006331", "ENSPREG00000006362", "ENSPREG00000006387", "ENSPREG00000006501",
                         "ENSPREG00000019738", "ENSPREG00000019635", "ENSPREG00000020098", "ENSPREG00000020098")

# I included (starts after ENSPREG00000006501): 
# dnajc25, Brcc3, KCNV2, and nars1 to check for duplicates 
combo_gene_list <- do.call(c, list(gene_list_autosome, gene_list_sexchromo))


# Check autosomal snps and their read depth:

gene_MF_depth <- read.table(file = "check_for_dupes/MF_snps_w_genes.txt", header = TRUE)

impact_gene_list <- lapply(combo_gene_list, function(p) {
  gene_MF_depth[grepl(p, gene_MF_depth$geneid), ]
})

names(impact_gene_list) <- combo_gene_list

# Pull out each gene to look at the M:F Read Depth and the SNPs:

movingaverage <- function (x, window) {
  ma <- roll_mean(x, window, fill = NA)
}

windowsize = 40

ENSPREG00000000149_df <- impact_gene_list[["ENSPREG00000000149"]]
ENSPREG00000000149_df$snp <- ifelse(ENSPREG00000000149_df$POS == 4040847, "impact", "normal")
ENSPREG00000000149_df <- ENSPREG00000000149_df[order(ENSPREG00000000149_df$snp, decreasing = T), ]

ENSPREG00000000149_perm <- subset(ENSPREG00000000149_df, snp != "impact")
mean(ENSPREG00000000149_perm$MF_ratio) # 1.108784
MFautopermute_149 <- replicate(1000,mean(sample(ENSPREG00000000149_perm$MF_ratio,windowsize,replace = T)))
MFautoI25cov_149 <- quantile(MFautopermute_149, c(.025, .5, .975))[[1]]
MFautoI975cov_149 <- quantile(MFautopermute_149, c(.025, .5, .975))[[3]]

ENSPREG00000000149_plot <- ggplot(ENSPREG00000000149_df, aes(x = POS)) + geom_point(aes(x = POS, y = MF_ratio, col = snp), size = 2) +
  scale_colour_manual(values = c("red", "grey70")) + theme_classic() + scale_y_continuous(lim = c(0, 2.0)) + 
  annotate("rect", xmin = 4036395, xmax = 4069745, ymax = MFautoI975cov_149, ymin = MFautoI25cov_149, fill="skyblue", alpha = 0.4) +
  geom_hline(aes(yintercept = 1.108784), col = "blue", linetype="dashed", size = 1) +
  theme(
    legend.position = 'right', strip.placement = 'outside', strip.background = eb, strip.text.y.left = element_text(angle = 0,size = 12),
    axis.text = element_text(size = 12), axis.title = element_text(size = 15), plot.title = element_text(hjust = 0.5, size = 15)
  ) + ylab("M:F Read Depth") +  xlab("Gene (bp)") + ggtitle("ENSPREG00000000149")
ENSPREG00000000149_plot

####

ENSPREG00000000425_df <- impact_gene_list[["ENSPREG00000000425"]]
ENSPREG00000000425_df$snp <- ifelse(ENSPREG00000000425_df$POS == 12036079 | ENSPREG00000000425_df$POS == 12036080, "impact", "normal")
ENSPREG00000000425_df <- ENSPREG00000000425_df[order(ENSPREG00000000425_df$snp, decreasing = T), ]

ENSPREG00000000425_perm <- subset(ENSPREG00000000425_df, snp != "impact")
mean(ENSPREG00000000425_perm$MF_ratio) # 1.13197
MFautopermute_425 <- replicate(1000,mean(sample(ENSPREG00000000425_perm$MF_ratio,windowsize,replace = T)))
MFautoI25cov_425 <- quantile(MFautopermute_425, c(.025, .5, .975))[[1]]
MFautoI975cov_425 <- quantile(MFautopermute_425, c(.025, .5, .975))[[3]]

ENSPREG00000000425_plot <- ggplot(ENSPREG00000000425_df, aes(x = POS)) + geom_point(aes(x = POS, y = MF_ratio, col = snp), size = 2) +
  scale_colour_manual(values = c("red", "grey70")) + theme_classic() + scale_y_continuous(lim = c(0, 2.0)) + 
  annotate("rect", xmin = 12035859, xmax = 12036488, ymax = MFautoI975cov_425, ymin = MFautoI25cov_425, fill="skyblue", alpha = 0.4) +
  geom_hline(aes(yintercept = 1.13197), col = "blue", linetype="dashed", size = 1) +
  theme(
    legend.position = 'right', strip.placement = 'outside', strip.background = eb, strip.text.y.left = element_text(angle = 0,size = 12),
    axis.text = element_text(size = 12), axis.title = element_text(size = 15), plot.title = element_text(hjust = 0.5, size = 15)
  ) + ylab("M:F Read Depth") +  xlab("Gene (bp)") + ggtitle("NANOS3")
ENSPREG00000000425_plot

####

ENSPREG00000000486_df <- impact_gene_list[["ENSPREG00000000486"]]
ENSPREG00000000486_df$snp <- ifelse(ENSPREG00000000486_df$POS == 13099938, "impact", "normal")
ENSPREG00000000486_df <- ENSPREG00000000486_df[order(ENSPREG00000000486_df$snp, decreasing = T), ]

ENSPREG00000000486_perm <- subset(ENSPREG00000000486_df, snp != "impact")
mean(ENSPREG00000000486_perm$MF_ratio) # 1.143596
MFautopermute_486 <- replicate(1000,mean(sample(ENSPREG00000000486_perm$MF_ratio,windowsize,replace = T)))
MFautoI25cov_486 <- quantile(MFautopermute_486, c(.025, .5, .975))[[1]]
MFautoI975cov_486 <- quantile(MFautopermute_486, c(.025, .5, .975))[[3]]

ENSPREG00000000486_plot <- ggplot(ENSPREG00000000486_df, aes(x = POS)) + geom_point(aes(x = POS, y = MF_ratio, col = snp), size = 2) +
  scale_colour_manual(values = c("red", "grey70")) + theme_classic() + scale_y_continuous(lim = c(0, 2.0)) +
  annotate("rect", xmin = 13061952, xmax = 13100864, ymax = MFautoI975cov_486, ymin = MFautoI25cov_486, fill="skyblue", alpha = 0.4) +
  geom_hline(aes(yintercept = 1.143596), col = "blue", linetype="dashed", size = 1) +
  theme(
    legend.position = 'right', strip.placement = 'outside', strip.background = eb, strip.text.y.left = element_text(angle = 0,size = 12),
    axis.text = element_text(size = 12), axis.title = element_text(size = 15), plot.title = element_text(hjust = 0.5, size = 15)
  ) + ylab("M:F Read Depth") +  xlab("Gene (bp)") + ggtitle("BCL11AB")
ENSPREG00000000486_plot

####

ENSPREG00000020299_df <- impact_gene_list[["ENSPREG00000020299"]]
ENSPREG00000020299_df$snp <- ifelse(ENSPREG00000020299_df$POS == 26058061, "impact", "normal")
ENSPREG00000020299_df <- ENSPREG00000020299_df[order(ENSPREG00000020299_df$snp, decreasing = T), ]

ENSPREG00000020299_perm <- subset(ENSPREG00000020299_df, snp != "impact")
mean(ENSPREG00000020299_perm$MF_ratio) # 0.9821289
MFautopermute_299 <- replicate(1000,mean(sample(ENSPREG00000020299_perm$MF_ratio,windowsize,replace = T)))
MFautoI25cov_299 <- quantile(MFautopermute_299, c(.025, .5, .975))[[1]]
MFautoI975cov_299 <- quantile(MFautopermute_299, c(.025, .5, .975))[[3]]

ENSPREG00000020299_plot <- ggplot(ENSPREG00000020299_df, aes(x = POS)) + geom_point(aes(x = POS, y = MF_ratio, col = snp), size = 2) +
  scale_colour_manual(values = c("red", "grey70")) + theme_classic() + scale_y_continuous(lim = c(0, 2.0)) +
  annotate("rect", xmin = 26052278, xmax = 26060042, ymax = MFautoI975cov_299, ymin = MFautoI25cov_299, fill="skyblue", alpha = 0.4) +
  geom_hline(aes(yintercept = 0.9821289), col = "blue", linetype="dashed", size = 1) +
  theme(
    legend.position = 'right', strip.placement = 'outside', strip.background = eb, strip.text.y.left = element_text(angle = 0,size = 12),
    axis.text = element_text(size = 12), axis.title = element_text(size = 15), plot.title = element_text(hjust = 0.5, size = 15)
  ) + ylab("M:F Read Depth") +  xlab("Gene (bp)") + ggtitle("SLC29A2")
ENSPREG00000020299_plot

#####

ENSPREG00000006517_df <- impact_gene_list[["ENSPREG00000006517"]]
ENSPREG00000006517_df$snp <- ifelse(ENSPREG00000006517_df$POS == 666979, "impact", "normal")
ENSPREG00000006517_df <- ENSPREG00000006517_df[order(ENSPREG00000006517_df$snp, decreasing = T), ]

ENSPREG00000006517_perm <- subset(ENSPREG00000006517_df, snp != "impact")
mean(ENSPREG00000006517_perm$MF_ratio) # 0.9718205
MFautopermute_6517 <- replicate(1000,mean(sample(ENSPREG00000006517_perm$MF_ratio,windowsize,replace = T)))
MFautoI25cov_6517 <- quantile(MFautopermute_6517, c(.025, .5, .975))[[1]]
MFautoI975cov_6517 <- quantile(MFautopermute_6517, c(.025, .5, .975))[[3]]

ENSPREG00000006517_plot <- ggplot(ENSPREG00000006517_df, aes(x = POS)) + geom_point(aes(x = POS, y = MF_ratio, col = snp), size = 2) +
  scale_colour_manual(values = c("red", "grey70")) + theme_classic() + scale_y_continuous(lim = c(0, 2.0)) +
  annotate("rect", xmin = 666593, xmax = 673446, ymax = MFautoI975cov_6517, ymin = MFautoI25cov_6517, fill="skyblue", alpha = 0.4) +
  geom_hline(aes(yintercept = 0.9718205), col = "blue", linetype="dashed", size = 1) +
  theme(
    legend.position = 'right', strip.placement = 'outside', strip.background = eb, strip.text.y.left = element_text(angle = 0,size = 12),
    axis.text = element_text(size = 12), axis.title = element_text(size = 15), plot.title = element_text(hjust = 0.5, size = 15)
  ) + ylab("M:F Read Depth") +  xlab("Gene (bp)") + ggtitle("TIMMDC1")
ENSPREG00000006517_plot

####

ENSPREG00000021017_df <- impact_gene_list[["ENSPREG00000021017"]]
ENSPREG00000021017_df$snp <- ifelse(ENSPREG00000021017_df$POS == 20150913 | ENSPREG00000021017_df$POS == 20150922, "impact", "normal")
ENSPREG00000021017_df <- ENSPREG00000021017_df[order(ENSPREG00000021017_df$snp, decreasing = T), ]

ENSPREG00000021017_perm <- subset(ENSPREG00000021017_df, snp != "impact")
mean(ENSPREG00000021017_perm$MF_ratio) # 1.020936
MFautopermute_21017 <- replicate(1000,mean(sample(ENSPREG00000021017_perm$MF_ratio,windowsize,replace = T)))
MFautoI25cov_21017 <- quantile(MFautopermute_21017, c(.025, .5, .975))[[1]]
MFautoI975cov_21017 <- quantile(MFautopermute_21017, c(.025, .5, .975))[[3]]

ENSPREG00000021017_plot <- ggplot(ENSPREG00000021017_df, aes(x = POS)) + geom_point(aes(x = POS, y = MF_ratio, col = snp), size = 2) +
  scale_colour_manual(values = c("red", "grey70")) + theme_classic() + scale_y_continuous(lim = c(0, 2.0)) +
  annotate("rect", xmin = 20132140, xmax = 20169144, ymax = MFautoI975cov_21017, ymin = MFautoI25cov_21017, fill="skyblue", alpha = 0.4) +
  geom_hline(aes(yintercept = 1.020936), col = "blue", linetype="dashed", size = 1) +
  theme(
    legend.position = 'right', strip.placement = 'outside', strip.background = eb, strip.text.y.left = element_text(angle = 0,size = 12),
    axis.text = element_text(size = 12), axis.title = element_text(size = 15), plot.title = element_text(hjust = 0.5, size = 15)
  ) + ylab("M:F Read Depth") +  xlab("Gene (bp)") + ggtitle("USP24")
ENSPREG00000021017_plot

####

ENSPREG00000001126_df <- impact_gene_list[["ENSPREG00000001126"]]
ENSPREG00000001126_df$snp <- ifelse(ENSPREG00000001126_df$POS == 3369883 | ENSPREG00000001126_df$POS == 3369885, "impact", "normal")
ENSPREG00000001126_df <- ENSPREG00000001126_df[order(ENSPREG00000001126_df$snp, decreasing = T), ]

ENSPREG00000001126_perm <- subset(ENSPREG00000001126_df, snp != "impact")
mean(ENSPREG00000001126_perm$MF_ratio) # 1.15532
MFautopermute_1126 <- replicate(1000,mean(sample(ENSPREG00000001126_perm$MF_ratio,windowsize,replace = T)))
MFautoI25cov_1126 <- quantile(MFautopermute_1126, c(.025, .5, .975))[[1]]
MFautoI975cov_1126 <- quantile(MFautopermute_1126, c(.025, .5, .975))[[3]]

ENSPREG00000001126_plot <- ggplot(ENSPREG00000001126_df, aes(x = POS)) + geom_point(aes(x = POS, y = MF_ratio, col = snp), size = 2) +
  scale_colour_manual(values = c("red", "grey70")) + theme_classic() + scale_y_continuous(lim = c(0, 2.0)) +
  annotate("rect", xmin = 3369212, xmax = 3386772, ymax = MFautoI975cov_1126, ymin = MFautoI25cov_1126, fill="skyblue", alpha = 0.4) +
  geom_hline(aes(yintercept = 1.15532), col = "blue", linetype="dashed", size = 1) +
  theme(
    legend.position = 'right', strip.placement = 'outside', strip.background = eb, strip.text.y.left = element_text(angle = 0,size = 12),
    axis.text = element_text(size = 12), axis.title = element_text(size = 15), plot.title = element_text(hjust = 0.5, size = 15)
  ) + ylab("M:F Read Depth") +  xlab("Gene (bp)") + ggtitle("ENSPREG00000001126")
ENSPREG00000001126_plot

#### 

ENSPREG00000012471_df <- impact_gene_list[["ENSPREG00000012471"]]
ENSPREG00000012471_df$snp <- ifelse(ENSPREG00000012471_df$POS == 5026983, "impact", "normal")
ENSPREG00000012471_df <- ENSPREG00000012471_df[order(ENSPREG00000012471_df$snp, decreasing = T), ]

ENSPREG00000012471_perm <- subset(ENSPREG00000012471_df, snp != "impact")
mean(ENSPREG00000012471_perm$MF_ratio) # 0.981656
MFautopermute_12471 <- replicate(1000,mean(sample(ENSPREG00000012471_perm$MF_ratio,windowsize,replace = T)))
MFautoI25cov_12471 <- quantile(MFautopermute_12471, c(.025, .5, .975))[[1]]
MFautoI975cov_12471 <- quantile(MFautopermute_12471, c(.025, .5, .975))[[3]]

ENSPREG00000012471_plot <- ggplot(ENSPREG00000012471_df, aes(x = POS)) + geom_point(aes(x = POS, y = MF_ratio, col = snp), size = 2) +
  scale_colour_manual(values = c("red", "grey70")) + theme_classic() + scale_y_continuous(lim = c(0, 2.0)) +
  annotate("rect", xmin = 5023901, xmax = 5037477, ymax = MFautoI975cov_12471, ymin = MFautoI25cov_12471, fill="skyblue", alpha = 0.4) +
  geom_hline(aes(yintercept = 0.981656), col = "blue", linetype="dashed", size = 1) +
  theme(
    legend.position = 'right', strip.placement = 'outside', strip.background = eb, strip.text.y.left = element_text(angle = 0,size = 12),
    axis.text = element_text(size = 12), axis.title = element_text(size = 15), plot.title = element_text(hjust = 0.5, size = 15)
  ) + ylab("M:F Read Depth") +  xlab("Gene (bp)") + ggtitle("ALDH16A1")
ENSPREG00000012471_plot

####

ENSPREG00000004135_df <- impact_gene_list[["ENSPREG00000004135"]]
ENSPREG00000004135_df$snp <- ifelse(ENSPREG00000004135_df$POS == 22642317, "impact", "normal")
ENSPREG00000004135_df <- ENSPREG00000004135_df[order(ENSPREG00000004135_df$snp, decreasing = T), ]

ENSPREG00000004135_perm <- subset(ENSPREG00000004135_df, snp != "impact")
mean(ENSPREG00000004135_perm$MF_ratio) # 1.126501
MFautopermute_4135 <- replicate(1000,mean(sample(ENSPREG00000004135_perm$MF_ratio,windowsize,replace = T)))
MFautoI25cov_4135 <- quantile(MFautopermute_4135, c(.025, .5, .975))[[1]]
MFautoI975cov_4135 <- quantile(MFautopermute_4135, c(.025, .5, .975))[[3]]

ENSPREG00000004135_plot <- ggplot(ENSPREG00000004135_df, aes(x = POS)) + geom_point(aes(x = POS, y = MF_ratio, col = snp), size = 2) +
  scale_colour_manual(values = c("red", "grey70")) + theme_classic() + scale_y_continuous(lim = c(0, 2.0)) +
  annotate("rect", xmin = 22639730, xmax = 22645801, ymax = MFautoI975cov_4135, ymin = MFautoI25cov_4135, fill="skyblue", alpha = 0.4) +
  geom_hline(aes(yintercept = 1.126501), col = "blue", linetype="dashed", size = 1) +
  theme(
    legend.position = 'right', strip.placement = 'outside', strip.background = eb, strip.text.y.left = element_text(angle = 0,size = 12),
    axis.text = element_text(size = 12), axis.title = element_text(size = 15), plot.title = element_text(hjust = 0.5, size = 15)
  ) + ylab("M:F Read Depth") +  xlab("Gene (bp)") + ggtitle("ENSPREG00000004135")
ENSPREG00000004135_plot

####

ENSPREG00000015247_df <- impact_gene_list[["ENSPREG00000015247"]]
ENSPREG00000015247_df$snp <- ifelse(ENSPREG00000015247_df$POS == 18673421, "impact", "normal")
ENSPREG00000015247_df <- ENSPREG00000015247_df[order(ENSPREG00000015247_df$snp, decreasing = T), ]

ENSPREG00000015247_perm <- subset(ENSPREG00000015247_df, snp != "impact")
mean(ENSPREG00000015247_perm$MF_ratio) # 1.063347
MFautopermute_15247 <- replicate(1000,mean(sample(ENSPREG00000015247_perm$MF_ratio,windowsize,replace = T)))
MFautoI25cov_15247 <- quantile(MFautopermute_15247, c(.025, .5, .975))[[1]]
MFautoI975cov_15247 <- quantile(MFautopermute_15247, c(.025, .5, .975))[[3]]

ENSPREG00000015247_plot <- ggplot(ENSPREG00000015247_df, aes(x = POS)) + geom_point(aes(x = POS, y = MF_ratio, col = snp), size = 2) +
  scale_colour_manual(values = c("red", "grey70")) + theme_classic() + scale_y_continuous(lim = c(0, 2.0)) +
  annotate("rect", xmin = 18663441, xmax = 18677212, ymax = MFautoI975cov_15247, ymin = MFautoI25cov_15247, fill="skyblue", alpha = 0.4) +
  geom_hline(aes(yintercept = 1.063347), col = "blue", linetype="dashed", size = 1) +
  theme(
    legend.position = 'right', strip.placement = 'outside', strip.background = eb, strip.text.y.left = element_text(angle = 0,size = 12),
    axis.text = element_text(size = 12), axis.title = element_text(size = 15), plot.title = element_text(hjust = 0.5, size = 15)
  ) + ylab("M:F Read Depth") +  xlab("Gene (bp)") + ggtitle("ENSPREG00000015247")
ENSPREG00000015247_plot

####

ENSPREG00000011757_df <- impact_gene_list[["ENSPREG00000011757"]]
ENSPREG00000011757_df$snp <- ifelse(ENSPREG00000011757_df$POS == 6298147, "impact", "normal")
ENSPREG00000011757_df <- ENSPREG00000011757_df[order(ENSPREG00000011757_df$snp, decreasing = T), ]

ENSPREG00000011757_perm <- subset(ENSPREG00000011757_df, snp != "impact")
mean(ENSPREG00000011757_perm$MF_ratio) # 1.032733
MFautopermute_11757 <- replicate(1000,mean(sample(ENSPREG00000011757_perm$MF_ratio,windowsize,replace = T)))
MFautoI25cov_11757 <- quantile(MFautopermute_11757, c(.025, .5, .975))[[1]]
MFautoI975cov_11757 <- quantile(MFautopermute_11757, c(.025, .5, .975))[[3]]
 
ENSPREG00000011757_plot <- ggplot(ENSPREG00000011757_df, aes(x = POS)) + geom_point(aes(x = POS, y = MF_ratio, col = snp), size = 2) +
  scale_colour_manual(values = c("red", "grey70")) + theme_classic() + scale_y_continuous(lim = c(0, 2.0)) +
  annotate("rect", xmin = 6297001, xmax = 6308472, ymax = MFautoI975cov_11757, ymin = MFautoI25cov_11757, fill="skyblue", alpha = 0.4) +
  geom_hline(aes(yintercept = 1.032733), col = "blue", linetype="dashed", size = 1) +
  theme(
    legend.position = 'right', strip.placement = 'outside', strip.background = eb, strip.text.y.left = element_text(angle = 0,size = 12),
    axis.text = element_text(size = 12), axis.title = element_text(size = 15), plot.title = element_text(hjust = 0.5, size = 15)
  ) + ylab("M:F Read Depth") +  xlab("Gene (bp)") + ggtitle("IFT74")
ENSPREG00000011757_plot

#####

ENSPREG00000013879_df <- impact_gene_list[["ENSPREG00000013879"]]
ENSPREG00000013879_df$snp <- ifelse(ENSPREG00000013879_df$POS == 15905017, "impact", "normal")
ENSPREG00000013879_df <- ENSPREG00000013879_df[order(ENSPREG00000013879_df$snp, decreasing = T), ]

ENSPREG00000013879_perm <- subset(ENSPREG00000013879_df, snp != "impact")
mean(ENSPREG00000013879_perm$MF_ratio) # 1.079553
MFautopermute_13879 <- replicate(1000,mean(sample(ENSPREG00000013879_perm$MF_ratio,windowsize,replace = T)))
MFautoI25cov_13879 <- quantile(MFautopermute_13879, c(.025, .5, .975))[[1]]
MFautoI975cov_13879 <- quantile(MFautopermute_13879, c(.025, .5, .975))[[3]]

ENSPREG00000013879_plot <- ggplot(ENSPREG00000013879_df, aes(x = POS)) + geom_point(aes(x = POS, y = MF_ratio, col = snp), size = 2) +
  scale_colour_manual(values = c("red", "grey70")) + theme_classic() + scale_y_continuous(lim = c(0, 2.0)) +
  annotate("rect", xmin = 15903334, xmax = 15906444, ymax = MFautoI975cov_13879, ymin = MFautoI25cov_13879, fill="skyblue", alpha = 0.4) +
  geom_hline(aes(yintercept = 1.079553), col = "blue", linetype="dashed", size = 1) +
  theme(
    legend.position = 'right', strip.placement = 'outside', strip.background = eb, strip.text.y.left = element_text(angle = 0,size = 12),
    axis.text = element_text(size = 12), axis.title = element_text(size = 15), plot.title = element_text(hjust = 0.5, size = 15)
  ) + ylab("M:F Read Depth") +  xlab("Gene (bp)") + ggtitle("SLC25A29")
ENSPREG00000013879_plot

####

ENSPREG00000004959_df <- impact_gene_list[["ENSPREG00000004959"]]
ENSPREG00000004959_df$snp <- ifelse(ENSPREG00000004959_df$POS == 15882247, "impact", "normal")
ENSPREG00000004959_df <- ENSPREG00000004959_df[order(ENSPREG00000004959_df$snp, decreasing = T), ]

ENSPREG00000004959_perm <- subset(ENSPREG00000004959_df, snp != "impact")
mean(ENSPREG00000004959_perm$MF_ratio) # 1.008537
MFautopermute_4959 <- replicate(1000,mean(sample(ENSPREG00000004959_perm$MF_ratio,windowsize,replace = T)))
MFautoI25cov_4959 <- quantile(MFautopermute_4959, c(.025, .5, .975))[[1]]
MFautoI975cov_4959 <- quantile(MFautopermute_4959, c(.025, .5, .975))[[3]]

ENSPREG00000004959_plot <- ggplot(ENSPREG00000004959_df, aes(x = POS)) + geom_point(aes(x = POS, y = MF_ratio, col = snp), size = 2) +
  scale_colour_manual(values = c("red", "grey70")) + theme_classic() + scale_y_continuous(lim = c(0, 2.0)) +
  annotate("rect", xmin = 15870593, xmax = 15883420, ymax = MFautoI975cov_4959, ymin = MFautoI25cov_4959, fill="skyblue", alpha = 0.4) +
  geom_hline(aes(yintercept = 1.008537), col = "blue", linetype="dashed", size = 1) +
  theme(
    legend.position = 'right', strip.placement = 'outside', strip.background = eb, strip.text.y.left = element_text(angle = 0,size = 12),
    axis.text = element_text(size = 12), axis.title = element_text(size = 15), plot.title = element_text(hjust = 0.5, size = 15)
  ) + ylab("M:F Read Depth") +  xlab("Gene (bp)") + ggtitle("WDR43")
ENSPREG00000004959_plot


library(patchwork)

(ENSPREG00000000149_plot | ENSPREG00000000425_plot | ENSPREG00000000486_plot | ENSPREG00000006517_plot) / 
  (ENSPREG00000021017_plot | ENSPREG00000001126_plot | ENSPREG00000012471_plot | ENSPREG00000004135_plot) /
  (ENSPREG00000015247_plot | ENSPREG00000011757_plot | ENSPREG00000013879_plot | ENSPREG00000004959_plot) + 
  plot_layout(guides = "collect")



###############################################################
#autosome_genes <- gene_MF_depth %>% filter(!geneid %in% combo_gene_list)
#autosome_genes <- subset(autosome_genes, CHROM != "LG12")
#autosome_genes$chrom_type <- "autosome"

##############################

avg_gene_MF <- read.table(file = "check_for_dupes/correct_ratio_MFDepth.txt", header = TRUE)

avg_gene_MF <- avg_gene_MF %>%
  mutate(chrom_type = ifelse(CHROM == "LG12", "sexchromo", "autosome"))

boxplot(avg_gene_MF$M_F_ratio ~ avg_gene_MF$chrom_type)

autosome_genes <- subset(avg_gene_MF, CHROM != "LG12")
summary(autosome_genes$M_F_ratio)

#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#  0.009901 1.011136 1.080857 1.081389 1.133643 2.996672 

IQR(autosome_genes$M_F_ratio) # 0.1225071
# Upper whisker = 1.5 * IQR = 0.18376065
# Median + Upper whisker = 1.26461765

mean_pute_genes <- avg_gene_MF %>%
  filter(geneid %in% combo_gene_list)

#check male limited SNPs and if they have significantly different M_F_ratio compared to the autosomal average
LG12_male_limited <- read.table(file = "check_for_dupes/LG12_male_limited_genes.txt", header = T)
autosomal_male_snps <- read.table(file = "Autosomal_male_limited_SNPs.txt", header = T)

################################

gene_list <- c(autosomal_male_snps$gene_id, LG12_male_limited$ensembl_ID)

autosome_genes <- autosome_genes %>% 
  filter(!(geneid %in% gene_list))
mean(autosome_genes$M_F_ratio)
# 1.081317

add_chrom_column <- function(dfA, list) {
  dfA$gene_detail <- ifelse(dfA$geneid %in% list, "pute_dupe", "autosomal")
  return(dfA)
}

org_gene_MFDepth <- add_chrom_column(gene_MF_depth, gene_list)
genes_to_test <- org_gene_MFDepth %>%
  filter(geneid %in% gene_list)

autosome_genes_values <- autosome_genes$M_F_ratio
pute_dupe_values <- genes_to_test$M_F_ratio

t.test(autosome_genes_values, pute_dupe_values, alternative = "two.sided")

#Output:

data:  autosome_genes_values and pute_dupe_values
t = -2.8263, df = 64.214, p-value = 0.00627
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  -0.08939293 -0.01535638
sample estimates:
  mean of x mean of y 
1.081317  1.133691

# Any gene that has a M:F Read Depth > 1.1333691 is likely a partial dupe

write.table(genes_to_test, "check_for_dupes/organize_the_pute_dupes.txt", quote = F, row.names = F, sep = "\t")

########


genes_to_test <- gene_MF_depth %>%
  filter(geneid %in% gene_list)

median(autosome_genes$M_F_ratio)
# 1.080857
summary(autosome_genes)
#Output:
#M:F ratio: min = 0.009901; 1st = 1.011136; median = 1.080857; mean = 1.081389; 3rd = 1.133643; max = 2.996672      
#upper whisker = min(max(x), Q_3 + 1.5 * IQR)
#lower whisker = max(min(x), Q_1 – 1.5 * IQR)

#Values above the upper whisker of autosomes = IQR = 0.122507 * 1.5 + Q3 = 1.3174035 



sexchromo_genes <- subset(gene_MF_depth, CHROM == "LG12")
median(sexchromo_genes$M_F_ratio)
# 1.094684
summary(sexchromo_genes)
#Output:
#M:F ratio: min = 0.822; 1st = 1.025; median = 1.095; mean = 1.094; 3rd = 1.137; max = 2.111
#upper whisker = min(max(x), Q_3 + 1.5 * IQR)
#lower whisker = max(min(x), Q_1 – 1.5 * IQR)


ggplot() +
  geom_density(aes(M_F_ratio, fill = "Autosome"), alpha = 0.4, data = autosome_genes) +
  geom_density(aes(M_F_ratio, fill = "Sex_chromosome"), alpha = 0.5, data = sexchromo_genes) +
  scale_fill_manual(name="Chromosome Type", 
                    values = c(Autosome = "grey", Sex_chromosome = "purple")) +
  theme_classic()


####### Everything below was not used to calculate gene duplicates #######

########### 2. Test for noramlity: 

ggqqplot(gene_MF_depth$M_F_ratio)
ggqqplot(sexchromo_genes$M_F_ratio)
ggqqplot(autosome_genes$M_F_ratio)





