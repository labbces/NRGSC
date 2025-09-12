wd = "/Storage/data1/jorge.munoz/NRGSC.new/networks/code/new"
setwd(wd)
# set working directory
#setwd(wd)
## load libraries 
library(DESeq2, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(dynamicTreeCut, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(fastcluster, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(WGCNA, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(reshape2, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
#library(DESeq2)
#library(WGCNA)
library(corpcor, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(longitudinal, lib.loc="/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(fdrtool, lib.loc="/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(GeneNet, lib.loc="/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(dplyr, lib.loc="/Storage/data1/jorge.munoz/NRGSC.new/libraries")
#library(dplyr)
#library(MCL)
library(MCL, lib.loc="/Storage/data1/jorge.munoz/NRGSC.new/libraries")
#library(igraph)
library(igraph, lib.loc="/Storage/data1/jorge.munoz/NRGSC.new/libraries")
#library(kernlab)
library(kernlab, lib.loc="/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(HiClimR, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries/")
#allowWGCNAThreads(5)
## variables
raw_results_path <- "/Storage/data1/jorge.munoz/NRGSC.new/networks/results/diffexpr_logfold_1.csv"
network_results_path <- "/Storage/data1/jorge.munoz/NRGSC.new/networks/results"
metadata <- read.table("/Storage/data1/jorge.munoz/NRGSC.new/data/metadata_complete.csv", sep=",", header = T)

## load raw data and aggregate for replicates
raw <- read.table(raw_results_path, header = T, sep = ",")
 filter genes with too many missing data
good <- goodSamplesGenes(raw, verbose = 3)
filtered <- as.matrix((raw[,good$goodGenes]))
# apply vst
vst <- varianceStabilizingTransformation(filtered)
# see distribution of vst
#hist <- hist(vst, xlim=c(2,20), breaks = 18)
#dev.off()
# filter genes with at least a value of expression > x among samples
# CORRER CORRELACIONES SIN FILTRAR POR VST
filtered_vst <- vst[apply(vst, 1, max) > 7,]
# calculate pearson correlation
#partial_cor <- ggm.estimate.pcor(t(filtered_vst))
#pcor <- cor(t(filtered_vst))
pcor <- fastCor(t(filtered_vst), nSplit = 40, upperTri = FALSE, optBLAS = TRUE, verbose = TRUE)
pcor[(abs(pcor) < 0.9)] <- 0
adj_list <- melt(pcor)
adj_filtered <- adj_list %>% filter(Var1 != Var2) %>% filter(value != 0)
write.table(adj_filtered, file = "./../../results/cool_adjacency_list.triples", row.names = F, col.names = F , quote = F, sep = " ")
#######################

# MCL in external sofware ./../../results/cool_adjacency_list.tsv

#######################


