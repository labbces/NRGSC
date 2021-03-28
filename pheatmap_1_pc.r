# script for ploting pretty heatmaps for contrast RNAseq data
# jmmunozp@usp.br
#install.packages("tidyverse")
#install.packages(c("pheatmap", "BiocManager"))
#BiocManager::install("DESeq2")

# variables
wd = "/home/j/BIOINFORMATICA/NRGSC_old"
mdata = "/metadata_complete.csv"
vst_path = "./DEA/complete/vst.Rdata"

# libraries
library(tidyverse)
library(DESeq2)
library(pheatmap)

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
pheatmap(dat, main =  paste( title," upregulated", sep = ""), annotation_col = anot, show_rownames = T, fontsize_row = 8, cellwidth = 20)
dev.off()

# downregulated genes
downregulated <- dea_contrast[head(order(dea_contrast$log2FoldChange),20),]
dat <- as.data.frame(t(contrast[downregulated$names]))
png(paste("./DEA/downregulated_", title, ".png", sep = ""),  width = 15*1.3, height = 15, res = 320, units = "cm")
pheatmap(dat, main = paste( title," downregulated", sep = ""), annotation_col = anot, show_rownames = T, fontsize_row = 8, cellwidth = 20)
dev.off()
}

# call function with different contrast
contrast_heatmaps(vst,c(34,30,26,46,42,38), "NR_270_B vs NR_10_B", "/home/j/BIOINFORMATICA/NRGSC_old/DEA/complete/filter_2NR_270_B_vs_NR_10_B.csv")




