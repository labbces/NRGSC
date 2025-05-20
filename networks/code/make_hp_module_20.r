library("igraph", lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
table <- read.table("/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/pearson_current.triples", header = F, sep = " ")
module20 <- read.table("/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/pearson_mcl/out.1.8.formated.csv", head = F, sep=" ", )

anot_row <- read.table("/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/TF/new/anot_row", sep = ",", header = T, colClasses = c("factor", "factor", "factor" , "factor", "numeric", "numeric"), row.names = 1)

module20_id <- module20[module20$V2==20,][1]

#colnames(table) <- c("V1", "V2", "weight")
#table$weight <- abs(table$weight)
#G <- simplify(graph_from_data_frame(table, directed = F), remove.multiple = T)

#degree <- degree(G)
#degree <- degree[module20_id$V1]
#betweenness <- betweenness(G))
#betweenness <- betweenness[module20_id$V1]

#family <- anot_row[module20_id,]


# Read nitrogen modules
#Nmods <- read.table("/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perGENOTYPE/pearson_GO_1.5/Nmod_full", header = F)
#colnames(Nmods) <- "Mod No"
tf <- read.table("/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/TF/new/TFs_in_network.ids")
#colnames(tf) <- "TF"
# files
#modules_path <- "./../../results/modules_formated.csv"
#modules_path <- "/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perGENOTYPE/pearson_mcl/out.1.8.formated.csv"
vst_path <- "/Storage/data1/jorge.munoz/NRGSC.new/networks/data/all_vst_counts.tsv"
metadata_path <-"/Storage/data1/jorge.munoz/NRGSC.new/networks/data/metadata_complete.csv"


#libraries
library(scico, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries/")
library(tidyverse, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries/")
library(DESeq2, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries/")
library(pheatmap, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries/")
library(viridis, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries/")
library(wesanderson, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries/")
library(RColorBrewer, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries/")

vst <- read.table(vst_path)
metadata <- read.table(metadata_path, sep = ",", header = T)

sample_table <- metadata
sample_table$Condition <- as.factor((sample_table$Condition))
sample_table$Genotype <- as.factor((sample_table$Genotype))
sample_table$DevStage <- as.factor((sample_table$DevStage))
sample_table$Individual <- as.factor((sample_table$Individual))
sample_table$Group <- as.factor(paste(sample_table$Genotype, '_', sample_table$Condition, '_', sample_table$DevStage, sep=''))

#df <- as.data.frame((assay(vst)))
annotation_col <- sample_table
anot <- select(annotation_col, "Group","DevStage","Genotype", "Condition")

df <- vst[module20_id$V1,]
df <- vst[tf$V1,]
rownames(anot) <- colnames(df)

# heat map with mean values for column i.e media por modulo
png("/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/TF_HP_samples_row.png", res = 300, width = 50, height = 40, units = "cm")
pheatmap(df,
         main ="Nitrogen Modules" ,
         scale = "row",
         annotation_col = anot,
         show_rownames = F,
         col = scico(100, palette = 'roma'),
         cluster_cols = T,
         cluster_rows = T,
         cellwidth = 12,
         cellheight = 12)
dev.off()

df2 <- data.frame(matrix(ncol = 16, nrow = dim(df)[1]))
colnames(df2) <- unique(anot$Group)
rownames(df2) <- rownames(df)

df2$NR_10_B <- (df$Sample_38 + df$Sample_42 + df$Sample_46)/3
df2$NR_10_B0 <- (df$Sample_37 + df$Sample_41 + df$Sample_45)/3
df2$NR_10_M <- (df$Sample_39 + df$Sample_43 + df$Sample_47)/3
df2$NR_10_P <- (df$Sample_40 + df$Sample_44 + df$Sample_48)/3
df2$NR_270_B <- (df$Sample_26 + df$Sample_30 + df$Sample_34)/3
df2$NR_270_B0 <- (df$Sample_25 + df$Sample_29 + df$Sample_33)/3
df2$NR_270_M <- (df$Sample_27 + df$Sample_31 + df$Sample_35)/3
df2$NR_270_P <- (df$Sample_28 + df$Sample_32 + df$Sample_36)/3
df2$R_10_B <- (df$Sample_14 + df$Sample_18 + df$Sample_22)/3
df2$R_10_B0 <- (df$Sample_13 + df$Sample_17 + df$Sample_21)/3
df2$R_10_M  <- (df$Sample_15 + df$Sample_19 + df$Sample_23)/3
df2$R_10_P <- (df$Sample_16 + df$Sample_20 + df$Sample_24)/3
df2$R_270_B <- (df$Sample_2 + df$Sample_6 + df$Sample_10)/3
df2$R_270_B0 <- (df$Sample_1 + df$Sample_5 + df$Sample_9)/3
df2$R_270_M <- (df$Sample_3 + df$Sample_7 + df$Sample_11)/3
df2$R_270_P <- (df$Sample_4 + df$Sample_8+ df$Sample_12)/3

write.table(df2, "/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/interesting_tf_matrix_for_pca.csv", sep = ",")
anot2 <- read.table("/Storage/data1/jorge.munoz/NRGSC.new/networks/data/metadata_complete_summarized.csv", sep = ",", header = T, colClasses = c("factor", "factor", "factor" ), row.names = 1)

#anot_row <- read.table("/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/TF/new/anot_row", sep = ",", header = T, colClasses = c("factor", "factor", "factor" , "factor", "numeric", "numeric"), row.names = 1)

anot_row <- read.table("/Storage/data1/jorge.munoz/NRGSC.new/networks/code/new/row_anot_mod20.csv", sep= ",", header = T, colClasses = c("factor", "factor", "factor", "factor", "numeric", "numeric"), row.names = 1)




#anot_row$Betweenness <- as.numeric(anot_row$Betweenness)
#anot_row$Degree <- as.numeric(anot_row$Degree)


colorDev=wes_palette("Zissou1", 4, type = "discrete")
colorGen=wes_palette("Darjeeling1", 2, type = "discrete")
colorCond=wes_palette("Chevalier1", 2, type = "discrete")
colorDeg=wes_palette("Zissou1", 100, type = "continuous")
colorBet=wes_palette("Zissou1", 100, type = "continuous")


#colorFamily=c("#bef7ff", "#c3ebed", "#c8dfda", "#cdd2c7", "#d2c6b3", "#d8b89e", "#ddaa89", "#e49b71", "#eb8a57", "#f37638", "#ff580a")
colorFamily=RColorBrewer::brewer.pal(n = 11, name = "Paired")
colorClass=c("#446455", "#FDD262", "black")
#colorModule=c("#8af73b", "#97ef3c", "#a3e73c", "#adde3c", "#b6d63d", "#becd3d", "#c5c53d", "#ccbc3e", "#d2b33e", "#d8a93e", "#dda03e", "#e2963e", "#e78b3e", "#eb803e", "#ef753e", "#f3683e", "#f65a3e", "#f9493d", "#fc343d", "#ff083d")
colorModule <- c("black","forestgreen", "red2", "orange", "cornflowerblue","#446455", "darkolivegreen4", "indianred1", "tan4", "darkblue", "mediumorchid1","firebrick4",  "yellowgreen", "lightsalmon", "tan3", "tan1", "darkgray", "wheat4", "#DDAD4B", "chartreuse")

ann_colors = list(
  Degree=colorDeg,
  Betweenness=colorBet,
  DevStage = c(B=colorDev[1], B0=colorDev[2], M=colorDev[3], P=colorDev[4]),
  Genotype = c(R = colorGen[1], NR = colorGen[2]),
  Condition = c("10" = colorCond[1], "270" = colorCond[2]),
  Family=c(MYB_related=colorFamily[10],
           MYB=colorFamily[11],
	   ARF=colorFamily[1],
	   AUX_IAA=colorFamily[2],
	   bHLH=colorFamily[3],
           bZIP=colorFamily[4],
           CCAAT=colorFamily[5],
           FHA=colorFamily[6],
           PHD=colorFamily[7],
           TRAF=colorFamily[8],
           WRKY=colorFamily[9],
	   no_tf="black"),
Module=c(
M102=colorModule[1],
M111=colorModule[2],
M119=colorModule[3],
M129=colorModule[4],
M140=colorModule[5],
M15=colorModule[6],
M17=colorModule[7],
M185=colorModule[8],
M191=colorModule[9],
M194=colorModule[10],
M195=colorModule[11],
M2=colorModule[12],
M20=colorModule[13],
M45=colorModule[14],
M46=colorModule[15],
M48=colorModule[16],
M75=colorModule[17],
M8=colorModule[18],
M9=colorModule[19],
M96=colorModule[20]
),
Class=c(OTR=colorClass[1], TFF=colorClass[2], no_tf=colorClass[3]),
cor_Genotype=colorDeg,
cor_Condition=colorDeg
)

trait_genotype <-as.numeric(c("0", "0", "0", "0", "0", "0", "0", "0", "1", "1", "1", "1", "1", "1", "1", "1"))
trait_nitrogen_condition <-as.numeric(c("1", "1", "1", "1","0", "0", "0", "0","1", "1", "1", "1","0", "0", "0", "0"))
trait_leaf_section <- as.numeric(c("1", "2", "3", "4", "1", "2", "3", "4", "1", "2", "3", "4", "1", "2", "3", "4"))
trait_R_up_response <- as.numeric(c("0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "1", "1", "1", "1"))


cor_genotype <- abs(cor(t(df2), trait_genotype, method = "spearman"))
colnames(cor_genotype) <- "cor_Genotype"
anot_row3  <- merge(anot_row, cor_genotype, by = 0)
rownames(anot_row3) <- anot_row3$Row.names
anot_row3 <- anot_row3[-1]


cor_condition <- abs(cor(t(df2), trait_nitrogen_condition, method = "spearman"))
colnames(cor_condition) <- "cor_Condition"
anot_row3  <- merge(anot_row3, cor_condition, by = 0)
rownames(anot_row3) <- anot_row3$Row.names
anot_row3 <- anot_row3[-1]

#cor_segment <- cor(t(df2), trait_leaf_section, method = "spearman")  
#colnames(cor_segment) <- "Leaf_segment"
#anot_row3  <- merge(anot_row3, cor_segment, by = 0)
#rownames(anot_row3) <- anot_row3$Row.names
#anot_row3 <- anot_row3[-1]

png("/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/module20_pretty_colors.png", width = 25, height = 35, units = "cm", res = 300)
pheatmap(df2,
         main ="",
         scale = "row",
         show_rownames = F,
         col = rev(scico(100, palette = 'roma')),
         cluster_cols = T,
         cluster_rows = T,
         annotation_col  = anot2,
         annotation_row = anot_row3,
         annotation_colors = ann_colors,
         annotation_names_row = T,
         cellwidth = 8,
         cellheight = 12
)
dev.off()
