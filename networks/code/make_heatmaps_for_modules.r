# Read nitrogen modules
Nmods <- read.table("/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/pearson_GO_1.8/all_mods.txt", header = F)
colnames(Nmods) <- "Mod No"

# files
#modules_path <- "./../../results/modules_formated.csv"
modules_path <- "/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/pearson_mcl/out.1.8.formated.csv"
vst_path <- "/Storage/data1/jorge.munoz/NRGSC.new/networks/data/all_vst_counts.tsv"
metadata_path <-"/Storage/data1/jorge.munoz/NRGSC.new/networks/data/metadata_complete.csv"


#libraries
library(scico, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries/")
library(tidyverse, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries/")
library(DESeq2, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries/")
library(pheatmap, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries/")
library(viridis, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries/")
library(wesanderson, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries/")

#setwd("/home/jmmunozp/things/heatmap_NRGSC/")
modules <- read.table(modules_path, row.names = 1, header = F)
colnames(modules) <- c("module_No")
# read full vst matrix 
vst <- read.table(vst_path)
metadata <- read.table(metadata_path, sep = ",", header = T)

# Make vectors for each intersting module
for (i in Nmods$'Mod No'){
assign(paste0("Module", i), modules %>% filter(module_No == i))
}

sample_table <- metadata
sample_table$Condition <- as.factor((sample_table$Condition))
sample_table$Genotype <- as.factor((sample_table$Genotype))
sample_table$DevStage <- as.factor((sample_table$DevStage))
sample_table$Individual <- as.factor((sample_table$Individual))
sample_table$Group <- as.factor(paste(sample_table$Genotype, '_', sample_table$Condition, '_', sample_table$DevStage, sep=''))

#df <- as.data.frame((assay(vst)))
annotation_col <- sample_table
anot <- select(annotation_col, "Group","DevStage","Genotype", "Condition")

for (i in Nmods$'Mod No'){
assign(paste0("module", i, "_dat"), vst[rownames(eval(as.name(paste0("Module",i)))),])
}

for (i in Nmods$'Mod No'){
assign(paste0("dat", i),as.data.frame(t(colMeans(eval(as.name(paste0("module",i, "_dat"))))))) 
labelling <-get(paste0("dat", i, sep = ""))
labels <- paste0("Module ", i, sep = "")
rownames(labelling) <- labels
assign(paste0("dat", i, sep= ""), labelling)
}

dat_list <- list()
for (i in Nmods$'Mod No'){
dat_list[[i]] <- eval(as.name(paste0("dat",i)))
}

df <- do.call("rbind", dat_list)
rownames(anot) <- colnames(df)
# heat map with mean values for column i.e media por modulo
png("/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/all_hp_color_ok_samples.png", res = 300, width = 7*800, height = 3*2850)
pheatmap(df,
         main ="Nitrogen Modules" ,
         scale = "row",
         annotation_col = anot,
         show_rownames = T,
         col = scico(100, palette = 'roma'),
         cluster_cols = T,
         cluster_rows = T,
         cellwidth = NA,
         cellheight = 8) 
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
df2$R_270_P <- (df$Sample_4 + df$Sample_8 + df$Sample_12)/3

#write.table(df2, "/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/all_matrix_for_pca.csv", sep = ",")
anot2 <- read.table("/Storage/data1/jorge.munoz/NRGSC.new/networks/data/metadata_complete_summarized.csv", sep = ",", header = T, colClasses = c("factor", "factor", "factor" ), row.names = 1)

colorDev=wes_palette("Zissou1", 4, type = "discrete")
colorGen=wes_palette("Darjeeling1", 2, type = "discrete")
colorCond=wes_palette("Chevalier1", 2, type = "discrete")

ann_colors = list(
  DevStage = c(B=colorDev[1], B0=colorDev[2], M=colorDev[3], P=colorDev[4]),
  Genotype = c(R = colorGen[1], NR = colorGen[2]),
  Condition = c("10" = colorCond[1], "270" = colorCond[2]))


df3 <- df2[1:30,]
png("/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/mean_modules_1_30.png", res = 300, width = 1500, height = 2850)
pheatmap(df3, 
         main ="Nitrogen Modules",
	 cellheight = 8,
	 cellwidth =NA,
         scale = "row",
         show_rownames = T,
         col = rev(scico(100, palette = 'roma')),
         cluster_cols = T,
         cluster_rows = T,
         annotation_col  = anot2,
         annotation_colors = ann_colors) 
dev.off()
