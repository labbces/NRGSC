load_cluster_libraries <- function(){
library(viridis, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(svglite, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(rlang, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(graph, lib.loc = "/Storage/data1/jorge.munoz/TOPGO/data")
library(GO.db, lib.loc = "/Storage/data1/jorge.munoz/TOPGO/data")
library(SparseM, lib.loc = "/Storage/data1/jorge.munoz/TOPGO/data")
library(topGO, lib.loc = "/Storage/data1/jorge.munoz/TOPGO/data")
library(dplyr, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(Rgraphviz, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(ggplot2, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(scales, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(clusterProfiler, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(parallel, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
}
setwd("/Storage/data1/jorge.munoz/NRGSC.new/networks/code/new")
load_cluster_libraries()
results_path_cluster <-"/Storage/data1/jorge.munoz/NRGSC.new/results/gene/GO_DEG/"
#dir.create(results_path_cluster)

#GO_ID = "/Storage/data1/jorge.munoz/NRGSC.new/networks/data/go_annotation_gene_ok.txt"
#UNIVERSE = "/Storage/data1/jorge.munoz/NRGSC.new/networks/data/go_anot.ids"
#INTERESTING = "/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/pearson_mcl/net.ids"
#OUT = "DEG"
GO_ID <- "/Storage/data1/jorge.munoz/NRGSC.new/networks/data/go_annotation_gene_ok.txt"
geneID2GO <- readMappings(file = GO_ID)

GO_ORA <-function(UNIVERSE, INTERESTING, OUT){
## read list of genes for all the network 
#geneID2GO <- readMappings(file = GO_ID)

geneNames <- read.table(UNIVERSE)
colnames(geneNames) <-"gene"

myInterestingGenes <- read.table(INTERESTING)
colnames(myInterestingGenes) <- "gene"

geneList <- as.factor(as.integer(geneNames$gene %in% myInterestingGenes$gene))
names(geneList) <- geneNames$gene

GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = geneList,
              annot = annFUN.gene2GO,
              gene2GO = geneID2GO)
allGO=usedGO(GOdata)

Classic <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
table <- GenTable(GOdata, Classic = Classic, topNodes = length(allGO), orderBy = 'Classic')

p.adj <- p.adjust(table$Classic,method="bonferroni")
tmp <- cbind(table,p.adj)
all_res_final <- tmp[which(tmp$p.adj<=0.001),]

#write.table(all_res_final, file = paste0(results_path_cluster, OUT, "_GO_TABLE", ".txt"), quote=FALSE, row.names=FALSE, sep = "\t")

ggdata <- all_res_final
aux <- go2term(all_res_final$GO.ID)
colnames(aux) <- c("GO.ID", "Lterm")

ggdata <- merge(ggdata, aux, by = "GO.ID")

ggdata$p.adj <- as.numeric(ggdata$p.adj)

ggdata <- ggdata[order(ggdata$p.adj),]
ggdata$Lterm <- factor(ggdata$Lterm, levels = rev(ggdata$Lterm)) # fixes order

write.table(all_res_final, file = paste0(results_path_cluster, OUT, "_GO_TABLE", ".txt"), quote=FALSE, row.names=FALSE, sep = "\t")
n <- 10
gg1 <- ggplot(head(ggdata,n),
aes(x = Lterm, y = -log10(p.adj) ))+		 
  
  expand_limits(y = 1) +
  geom_point(size = 6, colour = "black") +
  scale_size(range = c(2.5,12.5)) +
  
  xlab('') + ylab('-log(p)') +
  labs(
    title = OUT)+
  
  theme_bw() +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 42, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 32, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 32, face = 'bold', vjust = 1),    
    axis.text.x = element_text(angle = 0, size = 24, face = 'bold', hjust = 1.10, colour = 'black'),
    axis.text.y = element_text(angle = 0, size = 30, face = 'bold', vjust = 0.5, colour = 'black'),
    axis.title = element_text(size = 12, face = 'bold'),
    axis.title.x = element_text(size = 12, face = 'bold'),
    axis.title.y = element_text(size = 12, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold")) + 
  coord_flip()
ggsave(paste0(results_path_cluster, OUT, "_GO_PLOT.png"), plot = gg1, device = "png", width = 35, height = 30, dpi = 450, limitsize = F, units = "cm")
ggsave(paste0(results_path_cluster, OUT, "_GO_PLOT.svg"), plot = gg1, device = "svg", width = 35, height = 30, limitsize = F, units = "cm")
}
# make ORA for network
#GO_ORA("/Storage/data1/jorge.munoz/NRGSC.new/networks/data/go_annotation_gene_ok.txt", "/Storage/data1/jorge.munoz/NRGSC.new/networks/data/go_anot.ids", "/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/pearson_mcl/net.ids", "NET")
# make ORA for DEGs
#GO_ORA("/Storage/data1/jorge.munoz/NRGSC.new/networks/data/go_annotation_gene_ok.txt", "/Storage/data1/jorge.munoz/NRGSC.new/networks/data/go_anot.ids", "/Storage/data1/jorge.munoz/NRGSC.new/results/gene/TABLES/ids_diff_exp.ids","DEG")
# make ORA for subnetwork
#GO_ORA("/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/pearson_mcl/net.ids", "/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/TF/new/subnetwork_tf_ok.ids","SUBnet")

## MAKE ORA FOR CONSTRASTS
#GO_ORA("/Storage/data1/jorge.munoz/NRGSC.new/networks/data/go_anot.ids", "/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/names_R_270_P_vs_R_10_P.txt","R_P")
#GO_ORA("/Storage/data1/jorge.munoz/NRGSC.new/networks/data/go_anot.ids", "/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/names_R_270_M_vs_R_10_M.txt","R_M")
#GO_ORA("/Storage/data1/jorge.munoz/NRGSC.new/networks/data/go_anot.ids", "/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/names_R_270_B_vs_R_10_B.txt","R_B")
#GO_ORA("/Storage/data1/jorge.munoz/NRGSC.new/networks/data/go_anot.ids", "/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/names_R_270_B0_vs_R_10_B0.txt","R_B0")
#GO_ORA("/Storage/data1/jorge.munoz/NRGSC.new/networks/data/go_anot.ids", "/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/names_NR_270_P_vs_NR_10_P.txt","NR_P")
#GO_ORA("/Storage/data1/jorge.munoz/NRGSC.new/networks/data/go_anot.ids", "/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/names_NR_270_M_vs_NR_10_M.txt","NR_M")
#GO_ORA("/Storage/data1/jorge.munoz/NRGSC.new/networks/data/go_anot.ids", "/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/names_NR_270_B_vs_NR_10_B.txt","NR_B")
#GO_ORA("/Storage/data1/jorge.munoz/NRGSC.new/networks/data/go_anot.ids", "/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/names_NR_270_B0_vs_NR_10_B0.txt","NR_B0")

# MAKE ORA FOR DEG AND GENOTYPE
GO_ORA("/Storage/data1/jorge.munoz/NRGSC.new/networks/data/go_anot.ids", "/Storage/data1/jorge.munoz/NRGSC.new/results/gene/R_up_regulated.txt","R_up")
GO_ORA("/Storage/data1/jorge.munoz/NRGSC.new/networks/data/go_anot.ids", "/Storage/data1/jorge.munoz/NRGSC.new/results/gene/NR_up_regulated.txt","NR_up")
GO_ORA("/Storage/data1/jorge.munoz/NRGSC.new/networks/data/go_anot.ids", "/Storage/data1/jorge.munoz/NRGSC.new/results/gene/R_down_regulated.txt","R_down")
GO_ORA("/Storage/data1/jorge.munoz/NRGSC.new/networks/data/go_anot.ids", "/Storage/data1/jorge.munoz/NRGSC.new/results/gene/NR_down_regulated.txt","NR_down")
## MAKE ORA FOR CONTRASTS DOWN AND UP
## UP
#GO_ORA("/Storage/data1/jorge.munoz/NRGSC.new/networks/data/go_anot.ids", "/Storage/data1/jorge.munoz/NRGSC.new/results/gene/geneR_270_P_vs_R_10_P_up","R_P_up")
#GO_ORA("/Storage/data1/jorge.munoz/NRGSC.new/networks/data/go_anot.ids", "/Storage/data1/jorge.munoz/NRGSC.new/results/gene/geneR_270_M_vs_R_10_M_up","R_M_up")
#GO_ORA("/Storage/data1/jorge.munoz/NRGSC.new/networks/data/go_anot.ids", "/Storage/data1/jorge.munoz/NRGSC.new/results/gene/geneR_270_B_vs_R_10_B_up","R_B_up")
#GO_ORA("/Storage/data1/jorge.munoz/NRGSC.new/networks/data/go_anot.ids", "/Storage/data1/jorge.munoz/NRGSC.new/results/gene/geneR_270_B0_vs_R_10_B0_up","R_B0_up")
#GO_ORA("/Storage/data1/jorge.munoz/NRGSC.new/networks/data/go_anot.ids", "/Storage/data1/jorge.munoz/NRGSC.new/results/gene/geneNR_270_P_vs_NR_10_P_up","NR_P_up")
#GO_ORA("/Storage/data1/jorge.munoz/NRGSC.new/networks/data/go_anot.ids", "/Storage/data1/jorge.munoz/NRGSC.new/results/gene/NR_270_M_vs_NR_10_M_up","NR_M_up")
#GO_ORA("/Storage/data1/jorge.munoz/NRGSC.new/networks/data/go_anot.ids", "/Storage/data1/jorge.munoz/NRGSC.new/results/gene/geneNR_270_B_vs_NR_10_B_up","NR_B_up")
#GO_ORA("/Storage/data1/jorge.munoz/NRGSC.new/networks/data/go_anot.ids", "/Storage/data1/jorge.munoz/NRGSC.new/results/gene/NR_270_B0_vs_NR_10_B0_GO_up","NR_B0_up")
## DOWN
#GO_ORA("/Storage/data1/jorge.munoz/NRGSC.new/networks/data/go_anot.ids", "/Storage/data1/jorge.munoz/NRGSC.new/results/gene/geneR_270_P_vs_R_10_P_down","R_P_down")
#GO_ORA("/Storage/data1/jorge.munoz/NRGSC.new/networks/data/go_anot.ids", "/Storage/data1/jorge.munoz/NRGSC.new/results/gene/geneR_270_M_vs_R_10_M_down","R_M_down")
#GO_ORA("/Storage/data1/jorge.munoz/NRGSC.new/networks/data/go_anot.ids", "/Storage/data1/jorge.munoz/NRGSC.new/results/gene/geneR_270_B_vs_R_10_B_down","R_B_down")
#GO_ORA("/Storage/data1/jorge.munoz/NRGSC.new/networks/data/go_anot.ids", "/Storage/data1/jorge.munoz/NRGSC.new/results/gene/geneR_270_B0_vs_R_10_B0_down","R_B0_down")
#GO_ORA("/Storage/data1/jorge.munoz/NRGSC.new/networks/data/go_anot.ids", "/Storage/data1/jorge.munoz/NRGSC.new/results/gene/geneNR_270_P_vs_NR_10_P_down","NR_P_down")
#GO_ORA("/Storage/data1/jorge.munoz/NRGSC.new/networks/data/go_anot.ids", "/Storage/data1/jorge.munoz/NRGSC.new/results/gene/NR_270_M_vs_NR_10_M_down","NR_M_down")
#GO_ORA("/Storage/data1/jorge.munoz/NRGSC.new/networks/data/go_anot.ids", "/Storage/data1/jorge.munoz/NRGSC.new/results/gene/geneNR_270_B_vs_NR_10_B_down","NR_B_down")
#GO_ORA("/Storage/data1/jorge.munoz/NRGSC.new/networks/data/go_anot.ids", "/Storage/data1/jorge.munoz/NRGSC.new/results/gene/geneNR_270_B0_vs_NR_10_B0_GO_down","NR_B0_down")
