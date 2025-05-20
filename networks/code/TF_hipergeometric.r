#BiocManager::install("clusterProfiler", lib = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(clusterProfiler, lib.loc="/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(ggplot2, lib.loc="/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(svglite, lib.loc="/Storage/data1/jorge.munoz/NRGSC.new/libraries")


TFs_NETWORK <- "/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/TF/new/TF_families_in_NETWORK"
INTERESTING <- "/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/TF/new/TF_families_in_NETWORK"

genes <- as.matrix(read.table(TFs_NETWORK)[1])
colnames(genes) <- "gene_id"

test <- read.table(UNIVERSE,header = T)

all <- as.data.frame(cbind(test$term, test$gene))
colnames(all) <- c("term", "gene")

tf_erich <- enricher(gene = genes, TERM2GENE = all,  minGSSize = 5, maxGSSize = 10000, pvalueCutoff = 0.01)
#
results <- as.data.frame(cbind(tf_enrich$ID,tf_enrich$p.adjust))
colnames(results) <- c("term", "p.adjust")


