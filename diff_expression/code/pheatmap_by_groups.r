# variables
wd = "/Storage/data1/jorge.munoz/NRGSC.new/code"
# set working directory
setwd(wd)

# libraries
library(tidyverse, lib.loc = "./../libraries")
library(DESeq2, lib.loc = "./../libraries")
library(pheatmap, lib.loc = "./../libraries")
library(viridis, lib.loc = "./../libraries")

# function to graph heatmaps from contrast data

heatmaps <- function (vst_path, title, contrast_table, metadata) { 

  # read vst file
  df <- read.table(vst_path, header = T, sep = "\t")
  # read contrast table
  dea_contrast <- read.table(contrast_table, header = T, sep = "\t")
  # read metadata
  sample_table <-read.table(metadata, sep = "\t", header = T)
  sample_table$Condition <- as.factor((sample_table$Condition))
  sample_table$Genotype <- as.factor((sample_table$Genotype))
  sample_table$DevStage <- as.factor((sample_table$DevStage))
  sample_table$Individual <- as.factor((sample_table$Individual))
  sample_table$Group <- as.factor(paste(sample_table$Genotype, '_', sample_table$Condition, '_', sample_table$DevStage, sep=''))
  
  #df <- as.data.frame((assay(vst)))
  annotation_col <- sample_table
  anot <- select(annotation_col, "Group","DevStage","Genotype", "Condition")

  dat <- df[dea_contrast$names,]
  rownames(anot) <- colnames(dat)
  
  #heat map complete with all genes
  png(paste("./../results/hm_", title, ".png", sep = ""),  width = 15*1.3, height = 15, res = 320, units = "cm")
  pheatmap(dat, main = title, annotation_col = anot, show_rownames = F , cellwidth = 20, col = inferno(299))
  dev.off()

  # upregulated genes
 # upregulated <- dea_contrast[tail(order(dea_contrast$log2FoldChange),20),]
 # dat <- as.data.frame(df[upregulated$names,])
 # png(paste("./DEA/by_groups/upregulated_", title, ".png", sep = ""),  width = 15*1.3, height = 15, res = 320, units = "cm")
 # pheatmap(dat, main =  paste( title," upregulated", sep = ""), annotation_col = anot, show_rownames = T, fontsize_row = 8, cellwidth = 20, col = inferno(299))
 # dev.off()

  # downregulated genes
  #downregulated <- dea_contrast[head(order(dea_contrast$log2FoldChange),20),]
  #dat <- as.data.frame((df[downregulated$names,]))
  #png(paste("./DEA/by_groups/downregulated_", title, ".png", sep = ""),  width = 15*1.3, height = 15, res = 320, units = "cm")
  #pheatmap(dat, main = paste( title," downregulated", sep = ""), annotation_col = anot, show_rownames = T, fontsize_row = 8, cellwidth = 20, col = inferno(299))
  #dev.off()
}

# call function with different contrast

# NR_270_B0_vs_NR_10_B0
heatmaps("./../results/NR_270_B0-NR_10_B0-normalized_vst_counts.tsv", "NR_270_B0-NR_10_B0",  "./../results/NR_270_B0-NR_10_B0-DEA.tsv", "./../results/NR_270_B0-NR_10_B0-metadata.tsv")

# NR_270_B vs NR_10_B
heatmaps("./../results/NR_270_B-NR_10_B-normalized_vst_counts.tsv", "NR_270_B-NR_10_B",  "./../results/NR_270_B-NR_10_B-DEA.tsv", "./../results/NR_270_B-NR_10_B-metadata.tsv")

# NR_270_M_vs_NR_10_M
heatmaps("./../results/NR_270_M-NR_10_M-normalized_vst_counts.tsv", "NR_270_M-NR_10_M",  "./../results/NR_270_M-NR_10_M-DEA.tsv", "./../results/NR_270_M-NR_10_M-metadata.tsv")

# NR_270_P_vs_NR_10_P
heatmaps("./../results/NR_270_P-NR_10_P-normalized_vst_counts.tsv", "NR_270_P-NR_10_P",  "./../results/NR_270_P-NR_10_P-DEA.tsv", "./../results/NR_270_P-NR_10_P-metadata.tsv")

# R_270_B0_vs_R_10_B0
heatmaps("./../results/R_270_B0-R_10_B0-normalized_vst_counts.tsv", "R_270_B0-R_10_B0",  "./../results/R_270_B0-R_10_B0-DEA.tsv", "./../results/R_270_B0-R_10_B0-metadata.tsv")

# R_270_B_vs_R_10_B
heatmaps("./../results/R_270_B-R_10_B-normalized_vst_counts.tsv", "R_270_B-R_10_B",  "./../results/R_270_B-R_10_B-DEA.tsv", "./../results/R_270_B-R_10_B-metadata.tsv")

# R_270_M_vs_R_10_M
heatmaps("./../results/R_270_M-R_10_M-normalized_vst_counts.tsv", "R_270_M-R_10_M",  "./../results/R_270_M-R_10_M-DEA.tsv", "./../results/R_270_M-R_10_M-metadata.tsv")

# R_270_P_vs_R_10_P
heatmaps("./../results/R_270_P-R_10_P-normalized_vst_counts.tsv", "R_270_P-R_10_P",  "./../results/R_270_P-R_10_P-DEA.tsv", "./../results/R_270_P-R_10_P-metadata.tsv")
