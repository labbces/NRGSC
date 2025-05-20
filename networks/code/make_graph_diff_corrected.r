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
#raw <- aggregate(t(raw), by=list(Sample=metadata$Group), FUN=mean)
#rownames(raw) <- raw[,1]
#raw <- raw[-1]
# we need integers to calculate vst
#raw <- t(raw)
# filter genes with too many missing data
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
# ponerle un print a vst filtrado para ver cn cuantas columnas queda
# calculate pearson correlation
# fastcor MIRAR
#partial_cor <- ggm.estimate.pcor(t(filtered_vst))
#pcor <- cor(t(filtered_vst))
pcor <- fastCor(t(filtered_vst), nSplit = 40, upperTri = FALSE, optBLAS = TRUE, verbose = TRUE)
#pcor[(abs(partial_cor) < 0.0005)] <- 0
pcor[(abs(pcor) < 0.9)] <- 0
adj_list <- melt(pcor)
adj_filtered <- adj_list %>% filter(Var1 != Var2) %>% filter(value != 0)
#adj_filtered <- adj_list %>% filter(value !=0) %>% filter(value !=1)
#adj_filtered <- subset(adj_list, value !=0 & value !=1)
write.table(adj_filtered, file = "./../../results/cool_adjacency_list.triples", row.names = F, col.names = F , quote = F, sep = " ")
#######################

# MCL en softare externo ./../../results/cool_adjacency_list.tsv

#######################

#G <- graph_from_adjacency_matrix(pcor,  mode = "upper", weighted = TRUE, diag = FALSE)
#edge <- get.edgelist(G)
#edge$pcor <- E(G)$weight
#edge <- as.data.frame(edge)

# plot partial_cor to get cut
# exportar partial_cor para afinar parametros histograma 
#hist(edge$partial_cor, xlim = c(-1, 1), breaks = 1000, freq = F)
#hist(edge$partial_cor, xlim = c(-0.001, 0.001), breaks = 1000, freq = F)
#hist(edge$partial_cor, xlim = c(-0.05, 0.05), breaks = 1000, freq = F)
#dev.off()

#cool_edge_list <- edge %>% filter(abs(partial_cor) > 0.0005 && abs(pcor) > 0.35 )
#G <- delete_edges(G, which(abs(E(G)$weight) > 0.35 && abs(E(G)$partial_cor) > 0.0005)) 
# get pairs with abs(partial_cor) =>0.01 and pcor > 0.3
# COMPROBAR YO CON YO CORRELACION = 1	
# CAMBIAR LA DIAGONAL 
#list_network <- partial_cor[abs(pcor$Freq) > 0.35,] %>% filter(abs(Freq) > 0.0005) %>%  filter(abs(Freq) != 1)
#list_network_pcor <- pcor[abs(partial_cor$Freq) > 0.0005,] %>% filter(abs(Freq) > 0.35) %>% filter(abs(Freq) != 1)
# get adjacency matrix from edge list 

#G <- graph.data.frame(list_network_pcor, directed=FALSE)
#G <- graph.data.frame(cool_edge_list, directed=FALSE)
#A <- as_adjacency_matrix(G,type="both",names=TRUE, sparse=T, attr=" pcor" )
#$E(G)$weight <- list_network_pcor$Freq

#$T <- if(A != 0){A-1}
#A <- as_adjacency_matrix(G,type="both",names=TRUE, sparse=F)
#A[A > 0] <- 1
#A[A < 0] <- 1
#dim(A)
#A[1:10,1:10]
# write adjacency matrix
#write.table(as.matrix(A), "./../../results/adjacency_diff_expr_pcor.txt")
#write list netwokr
#write.table(list_network, "./../../results/list.network_pcor.txt")

###################################
# get modules with markov chain clustering
# cluster mcl
#mc <- mcl(pcor, addLoops = F, allow1 = T)
#mc
# COLOCAR EN 
#mc$K
#mccluster <- as.data.frame(mc$Cluster, optional = F)
# idenfy modules with spectral clustering
#modules <- as.data.frame(specc(as.matrix(pcor), mc$K))
#modules <- as.data.frame(clusters$Cluster)
# COMPARAR MODULOS POR MCL Y SPECC
#colnames(mccluster) <- "module_No"
# write modules.txt spectral clustering
#rite.table(modules, file = "./../../results/modules_specc_diff_expr_pcor.txt")
# write modules markov chain clustering 
#rownames(mccluster) <- rownames(pcor)
#colnames(mccluster) <- "module_No"
#write.table(mccluster, file = "./../../results/modules_mcc_diff_expr_pcor.txt")

#a<- read.table("./../../results/modules.txt")
#MA <- read.table("./../../results/adjacency.txt")

