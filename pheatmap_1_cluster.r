# script for ploting pretty heatmaps for contrast RNAseq data
# jmmunozp@usp.br
#install.packages("tidyverse")
#install.packages("BiocManager")
install.packages("pheatmap",
lib = "./libraries",
repos = 'http://cran.us.r-project.org')
#BiocManager::install("DESeq2")

# variables
wd = "/Storage/data1/jorge.munoz/NRGSC"
mdata = "./metadata_complete.csv"
vst_path = "./DEA/vst.Rdata"

# libraries
library(tidyverse, lib.loc = "./libraries")
library(DESeq2, lib.loc = "./libraries")
library(pheatmap, lib.loc = "./libraries")

# set working directory
setwd(wd)

sample_table <-read.table(mdata, sep = ",", header = T)
sample_table$Condition <- as.factor((sample_table$Condition))
sample_table$Genotype <- as.factor((sample_table$Genotype))
sample_table$DevStage <- as.factor((sample_table$DevStage))
sample_table$Individual <- as.factor((sample_table$Individual))
sample_table$Group <- as.factor(paste(sample_table$Genotype, '_', sample_table$Condition, '_', sample_table$DevStage, sep=''))

# read vst file
load(vst_path)
# function to graph heatmaps from contrast data

contrast_heatmaps <- function (vst, contrast_samples_index, title, contrast_table) { 

dea_contrast <- read.table(contrast_table, header = T, sep = ",")

df <- as.data.frame(t(assay(vst)))
# NR_270_B (Sample_34, Sample_30, Sample_26)
# NR_10_B (Sample_46, Sample_42, Sample_38)
contrast <- df[contrast_samples_index,]

annotation_col <- sample_table[contrast_samples_index,]
anot <- select(annotation_col, "Group","DevStage","Genotype", "Condition")
dat <- as.data.frame(t(contrast[dea_contrast$names]))
rownames(anot) <- colnames(dat)
#heat map complete with all genes
png(paste("./DEA/hm_", title, ".png", sep = ""),  width = 15*1.3, height = 15, res = 320, units = "cm")
pheatmap(dat, main = title, annotation_col = anot, show_rownames = F , cellwidth = 20)
dev.off()

# upregulated genes
upregulated <- dea_contrast[tail(order(dea_contrast$log2FoldChange),20),]
dat <- as.data.frame(t(contrast[upregulated$names]))
png(paste("./DEA/upregulated_", title, ".png", sep = ""),  width = 15*1.3, height = 15, res = 320, units = "cm")
pheatmap(dat, main =  paste( title," upregulated", sep = ""), annotation_col = anot, show_rownames = T, fontsize_row = 8, cellwidth = 20,)
dev.off()

# downregulated genes
downregulated <- dea_contrast[head(order(dea_contrast$log2FoldChange),20),]
dat <- as.data.frame(t(contrast[downregulated$names]))
png(paste("./DEA/downregulated_", title, ".png", sep = ""),  width = 15*1.3, height = 15, res = 320, units = "cm")
pheatmap(dat, main = paste( title," downregulated", sep = ""), annotation_col = anot, show_rownames = T, fontsize_row = 8, cellwidth = 20)
dev.off()
}

# call function with different contrast
# NR_270_B vs NR_10_B
contrast_heatmaps(vst,c(26,30,34,38,42,46), "NR_270_B vs NR_10_B", "./DEA/filter3_NR_270_B_vs_NR_10_B.csv")
# NR_270_B0_vs_NR_10_B0
contrast_heatmaps(vst,c(25,29,33,37,41,45), "NR_270_B0 vs NR_10_B0", "./DEA/filter3_NR_270_B0_vs_NR_10_B0.csv")
# NR_270_M_vs_NR_10_M
contrast_heatmaps(vst,c(27,31,35,39,43,47), "NR_270_M vs NR_10_M", "./DEA/filter3_NR_270_M_vs_NR_10_M.csv")
# NR_270_P_vs_NR_10_P
contrast_heatmaps(vst,c(28,32,36,40,44,48), "NR_270_P vs NR_10_P", "./DEA/filter3_NR_270_P_vs_NR_10_P.csv")
# R_270_B0_vs_R_10_B0
contrast_heatmaps(vst,c(1,5,9,13,17,21), "R_270_B0 vs R_10_B0", "./DEA/filter3_R_270_B0_vs_R_10_B0.csv")
# R_270_B_vs_R_10_B
contrast_heatmaps(vst,c(2,6,10,14,18,22), "R_270_B vs R_10_B", "./DEA/filter3_R_270_B_vs_R_10_B.csv")
# R_270_M_vs_R_10_M
contrast_heatmaps(vst,c(3,7,11,15,19,23), "R_270_M vs R_10_M", "./DEA/filter3_R_270_M_vs_R_10_M.csv")
# R_270_P_vs_R_10_P
contrast_heatmaps(vst,c(4,8,12,16,20,24), "R_270_P vs R_10_P", "./DEA/filter3_R_270_P_vs_R_10_P.csv")
