# set working directory
setwd(".")
## install packages
# install.packages(c("BiocManager","readr","dplyr", "magrittr","ggplot2","backports","tidyverse"),
# lib = "./libraries",
# repos = 'http://cran.us.r-project.org')


# BiocManager::install("tximport", lib = "./libraries")
# BiocManager::install("DESeq2", lib = "./libraries")
# BiocManager::install("RUVSeq", lib = "./libraries")
# BiocManager::install("EDASeq", lib = "./libraries")
# BiocManager::install("edgeR", lib = "./libraries")


#  load libraries

library(readr, lib.loc =  "./libraries")
library(dplyr, lib.loc =  "./libraries/")
library(magrittr,  lib.loc = "./libraries")
library(tximport, lib.loc = "./libraries")
library(DESeq2, lib.loc = "./libraries")
library(ggplot2, lib.loc =  "./libraries")
library(EDASeq, lib.loc =  "./libraries")
library(RUVSeq, lib.loc =  "./libraries")
library(edgeR, lib.loc = "./libraries")
library(backports, lib.loc = "./libraries")
library(tidyverse, lib.loc = "./libraries")

 
 #library(readr)
 #library(dplyr)
 #library(magrittr)
 #library(tximport)
 #library(DESeq2)
 #library(ggplot2)
 #library(EDASeq)
 #library(RUVSeq)
 #library(edgeR)
 #load metadata

sample_table <-read.table("metadata_complete.csv", sep = ",", header = T)
# delete_samples
#sample_table<-rows_delete(sample_table, tibble(Sample.Number = 37))
# load files paths
sample_files = paste0(pull(sample_table , "Sample_file"), "/quant.sf")
# name table columns
names(sample_files) = pull(sample_table, "Sample_file")
# relate genes to transcripts
tx2gene = read.table("txgen2.txt", sep = ",", col.names =c("genid","transid"))
# import count data to tximport
count_data = tximport( files = sample_files,
          type = "salmon",
          tx2gene =  tx2gene,
          ignoreTxVersion = T)

sample_table$Condition <- as.factor((sample_table$Condition))
sample_table$Genotype <- as.factor((sample_table$Genotype))
sample_table$DevStage <- as.factor((sample_table$DevStage))
sample_table$Individual <- as.factor((sample_table$Individual))

sample_table$Group <- as.factor(paste(sample_table$Genotype, '_', sample_table$Condition, '_', sample_table$DevStage, sep=''))

raw <- DESeqDataSetFromTximport(txi = count_data,
                                 colData = sample_table,
                                 design = ~ Group)

dim(raw)
keep<- rowSums(counts(raw))>0
table(keep)
raw <- raw[keep,]
dim(raw)

data <- estimateSizeFactors(raw)
data <- estimateDispersions(data) 
data <- nbinomWaldTest(data)

#### Differencial expression analyses
dea <- DESeq(data)
#resultsNames(data)

dea_analysis <- function( file , cond1, cond2, title, graphname ) {
        dea_contrast <- results(file, contrast=c("Group", cond1, cond2))
        dea_df <- as.data.frame(dea_contrast)
### filter DEGs 
        filter_1 <- dea_df[complete.cases(dea_df),]
        filter_2 <- filter_1[filter_1$padj < 0.05,]
        filter_3 <- filter_2[abs(filter_2$log2FoldChange) > 1,]
## volcano plot
        filter_1$test = filter_1$padj < 0.05 & abs(filter_1$log2FoldChange) > 1
        filter_1 = rownames_to_column(filter_1, var='ensgene')
        g = ggplot(filter_1, aes(x=log2FoldChange,
                                 y=-log10(padj), 
                                name=ensgene)) +
                geom_point(aes(colour=test), size=1, alpha=0.3) +
                scale_colour_manual(values=c('black', 'red')) +
                geom_vline(xintercept=1, colour='green', linetype=3) +
                geom_vline(xintercept=-1, colour='green', linetype=3) +
                geom_hline(yintercept=-log10(0.05), colour='blue', linetype=3) +
                labs(title = title) +
                theme_bw() +
                theme( text = element_text(size=22), 
                       panel.border = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.line = element_line(colour = "black"),
                       legend.position = 'none')
ggsave( g, filename = paste("./DEA/volcanoplot", title, ".png") ,units = "cm",width = 15*1.3, height = 15,dpi = 320)
##
filter_1$names <- row.names(filter_1)
filter_2$names <- row.names(filter_2)
filter_3$names <- row.names(filter_3)
#plotCounts(dea, gene = "NewTr2475572" ,intgroup = "Condition")
write_csv(filter_3, file = paste("./DEA/", "filter_3", title ,".csv", sep = ""))
write_csv(filter_2, file = paste("./DEA/", "filter_2", title ,".csv", sep = ""))
write_csv(filter_1, file = paste("./DEA/", "filter_1", title ,".csv", sep = ""))
}

dea_analysis(dea, "NR_270_B", "NR_10_B", "NR_270_B_vs_NR_10_B", "NR_270_B_vs_NR_10_B" )
dea_analysis(dea, "NR_270_B0", "NR_10_B0", "NR_270_B0_vs_NR_10_B0", "NR_270_B0_vs_NR_10_B0" )
dea_analysis(dea, "NR_270_M", "NR_10_M", "NR_270_M_vs_NR_10_M", "NR_270_M_vs_R_10_M" )
dea_analysis(dea, "NR_270_P", "NR_10_P", "NR_270_P_vs_NR_10_P", "NR_270_P_vs_R_10_P" )
dea_analysis(dea, "R_270_B0", "R_10_B0", "R_270_B0_vs_R_10_B0", "R_270_B0_vs_R_10_B0" )
dea_analysis(dea, "R_270_B", "R_10_B", "R_270_B_vs_R_10_B", "R_270_B_vs_R_10_B" )
dea_analysis(dea, "R_270_M", "R_10_M", "R_270_M_vs_R_10_M", "R_270_M_vs_R_10_M" )
dea_analysis(dea, "R_270_P", "R_10_P", "R_270_P_vs_R_10_P", "R_270_P_vs_R_10_P" )

save.image( file = "./DEA/DEA_volcanoplot.Rdata")
