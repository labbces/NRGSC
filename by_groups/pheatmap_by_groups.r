# variables
wd = "/Storage/data1/jorge.munoz/NRGSC"
# set working directory
setwd(wd)

# libraries
library(tidyverse, lib.loc = "./libraries")
library(DESeq2, lib.loc = "./libraries")
library(pheatmap, lib.loc = "./libraries")
library(viridis, lib.loc = "./libraries")


# function to graph heatmaps from contrast data

heatmaps <- function (vst_path, title, contrast_table) { 

  # read vst file
  load(vst_path)
  
  sample_table <-read.table("./metadata_complete.csv", sep = ",", header = T)
  sample_table$Condition <- as.factor((sample_table$Condition))
  sample_table$Genotype <- as.factor((sample_table$Genotype))
  sample_table$DevStage <- as.factor((sample_table$DevStage))
  sample_table$Individual <- as.factor((sample_table$Individual))
  sample_table$Group <- as.factor(paste(sample_table$Genotype, '_', sample_table$Condition, '_', sample_table$DevStage, sep=''))
  sample_table <- sample_table[vst$Sample.Number,]
  
  dea_contrast <- read.table(contrast_table, header = T, sep = ",")
  df <- as.data.frame((assay(vst)))
  annotation_col <- sample_table
  anot <- select(annotation_col, "Group","DevStage","Genotype", "Condition")

  dat <- df[dea_contrast$names,]

  rownames(anot) <- colnames(dat)
  #heat map complete with all genes
  png(paste("./DEA/by_groups/hm_", title, ".png", sep = ""),  width = 15*1.3, height = 15, res = 320, units = "cm")
  pheatmap(dat, main = title, annotation_col = anot, show_rownames = F , cellwidth = 20, col = inferno(299))
  dev.off()

  # upregulated genes
  upregulated <- dea_contrast[tail(order(dea_contrast$log2FoldChange),20),]
  dat <- as.data.frame(df[upregulated$names,])
  png(paste("./DEA/by_groups/upregulated_", title, ".png", sep = ""),  width = 15*1.3, height = 15, res = 320, units = "cm")
  pheatmap(dat, main =  paste( title," upregulated", sep = ""), annotation_col = anot, show_rownames = T, fontsize_row = 8, cellwidth = 20, col = inferno(299))
  dev.off()

  # downregulated genes
  downregulated <- dea_contrast[head(order(dea_contrast$log2FoldChange),20),]
  dat <- as.data.frame((df[downregulated$names,]))
  png(paste("./DEA/by_groups/downregulated_", title, ".png", sep = ""),  width = 15*1.3, height = 15, res = 320, units = "cm")
  pheatmap(dat, main = paste( title," downregulated", sep = ""), annotation_col = anot, show_rownames = T, fontsize_row = 8, cellwidth = 20, col = inferno(299))
  dev.off()
}

# call function with different contrast

# NR_270_B0_vs_NR_10_B0
heatmaps("./DEA/by_groups/NR_270_B0_vs_NR_10_B0_vst.Rdata", "NR_270_B0 vs NR_10_B0",  "DEA/by_groups/filter_3NR_270_B0_vs_NR_10_B0.csv")

# NR_270_B vs NR_10_B
heatmaps("./DEA/by_groups/NR_270_B_vs_NR_10_B_vst.Rdata", "NR_270_B vs NR_10_B",  "DEA/by_groups/filter_3NR_270_B_vs_NR_10_B.csv")

# NR_270_M_vs_NR_10_M
heatmaps("./DEA/by_groups/NR_270_M_vs_NR_10_M_vst.Rdata", "NR_270_M vs NR_10_M",  "DEA/by_groups/filter_3NR_270_M_vs_NR_10_M.csv")

# NR_270_P_vs_NR_10_P
heatmaps("./DEA/by_groups/NR_270_P_vs_NR_10_P_vst.Rdata", "NR_270_P vs NR_10_P",  "DEA/by_groups/filter_3NR_270_P_vs_NR_10_P.csv")

# R_270_B0_vs_R_10_B0
heatmaps("./DEA/by_groups/R_270_B0_vs_R_10_B0_vst.Rdata", "R_270_B0 vs R_10_B0",  "DEA/by_groups/filter_3R_270_B0_vs_R_10_B0.csv")

# R_270_B_vs_R_10_B
heatmaps("./DEA/by_groups/R_270_B_vs_R_10_B_vst.Rdata", "R_270_B vs R_10_B",  "DEA/by_groups/filter_3R_270_B_vs_R_10_B.csv")

# R_270_M_vs_R_10_M
heatmaps("./DEA/by_groups/R_270_M_vs_R_10_M_vst.Rdata", "R_270_M vs R_10_M",  "DEA/by_groups/filter_3R_270_M_vs_R_10_M.csv")

# R_270_P_vs_R_10_P
heatmaps("./DEA/by_groups/R_270_P_vs_R_10_P_vst.Rdata", "R_270_P vs R_10_P",  "DEA/by_groups/filter_3R_270_P_vs_R_10_P.csv")
