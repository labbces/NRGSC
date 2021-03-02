# set number of genes for empirical control for removing unwanted variation
EC = 100
# set working directory
setwd("/home/j/BIOINFORMATICA/NRGSC_old/")
## install packages
# install.packages(c("BiocManager","readr","dplyr", "magrittr","ggplot2"),
# lib = "./libraries",
# repos = 'http://cran.us.r-project.org')
# BiocManager::install("tximport", lib = "./libraries")
# BiocManager::install("DESeq2", lib = "./libraries")
# BiocManager::install("RUVSeq", lib = "./libraries")
# BiocManager::install("EDASeq", lib = "./libraries")
# BiocManager::install("edgeR", lib = "./libraries")
library(readr, lib.loc =  "./libraries")
library(dplyr, lib.loc =  "./libraries/")
library(magrittr,  lib.loc = "./libraries")
library(tximport, lib.loc = "./libraries")
library(DESeq2, lib.loc = "./libraries")
library(ggplot2, lib.loc =  "./libraries")
library(hexbin, lib.loc =  "./libraries")
library(EDASeq, lib.loc =  "./libraries")
library(RUVSeq, lib.loc =  "./libraries")
library(edgeR, lib.loc = "./libraries")
### function plot pca from deseq object
plot_pca_deseq <- function(data, shape, title, name, path) {
a <- plotPCA(data, intgroup ='Group') +
theme_bw() +
geom_point(size=4.5, aes( shape = shape))+
scale_shape_manual(values = seq(0,8)) +
labs(title = title, col="Group", shape = name) +
theme( text = element_text(size=22),
panel.border = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
plot.title = element_text(hjust = 0.5),
axis.line = element_line(colour = "black") )
ggsave( a, filename = path ,units = "cm",width = 15*1.3, height = 15,dpi = 320)
return(a)
}
plot_pca <- function(data, shape, title, percentage, path, legend, color) {
r <- ggplot(data, aes(x=PC1, y=PC2, colour=color, shape=shape))
r <- r + geom_point(size=4.5)
r <- r + theme_bw()
r <- r + labs(title = title, col="Group", shape= legend)
r <- r + xlab((paste0("PC1 : ",round(100 * (percentage[1]^2/sum(percentage^2)),2),"%")))
r <- r + ylab((paste0("PC2 : ",round(100 * (percentage[2]^2/sum(percentage^2)),2),"%")))
r <- r + theme( text = element_text(size=22),
panel.border = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_line(colour = "black"))
ggsave( r, filename = path ,units = "cm",width = 15*1.3, height = 15,dpi = 320)
return(r)
}
sample_table <-read.table("metadata_complete.csv", sep = ",", header = T)
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
sample_table$Group <- paste(sample_table$Genotype,'_',
sample_table$Condition,'_',
sample_table$DevStage, sep='')
sample_table$Group<-as.factor(sample_table$Group)
raw <- DESeqDataSetFromTximport(txi = count_data,
colData = sample_table,
design = ~ Group)
dim(raw)
keep<- rowSums(counts(raw))>0
table(keep)
raw <- raw[keep,]
dim(raw)
data <- estimateSizeFactors(raw)
vst <- varianceStabilizingTransformation(data)
## pca raw
raw_pca <-assay(data)
pca_raw <- prcomp(raw_pca,
center = T,
scale. = T)
pca_raw_data <- as.data.frame(pca_raw$rotation)
## plot raw by devstage
plot_pca(pca_raw_data, sample_table$DevStage, "RAW DATA BY DEVSTAGE", pca_raw$sdev,  "./raw_by_devstage.png", "DevStage", sample_table$Group)
## plot raw Genotype
plot_pca(pca_raw_data, sample_table$Genotype, "RAW DATA BY GENOTYPE", pca_raw$sdev,  "./raw_by_genotype.png", "Genotype", sample_table$Group)
## plot raw condition
plot_pca(pca_raw_data, sample_table$Condition, "RAW DATA BY CONDITION", pca_raw$sdev, "./raw_by_condition.png", "Condition", sample_table$Group)

## plot vst Devstage
plot_pca_deseq(vst, vst$DevStage, "VST BY DEVSTAGE", "Devstage",  "./vst/vst_devstage.png")
## plot vst Genotype
plot_pca_deseq(vst, vst$Genotype, "VST BY GENOTYPE", "Genotype",  "./vst/vst_genotype.png")
## plot vst Condition
plot_pca_deseq(vst, vst$Condition, "VST BY CONDITION", "Condition",  "./vst/vst_condition.png")

# vst R condition
vst_R <- vst[ ,vst$Genotype ==  "R" ]
plot_pca_deseq(vst_R, vst_R$Condition, "VST R BY CONDITION", "Condition",  "./vst/vst_R_condition.png")
# vst R devstage
plot_pca_deseq(vst_R, vst_R$DevStage, "VST R BY DEVSTAGE", "Devstage",  "./vst/vst_R_devstage.png")
# vst NR condition
vst_NR <- vst[ ,vst$Genotype ==  "NR" ]
plot_pca_deseq(vst_NR, vst_NR$Condition, "VST NR BY CONDITION", "Condition",  "./vst/vst_NR_condition.png")
# vst NR devstage
plot_pca_deseq(vst_NR, vst_R$DevStage, "VST NR BY DEVSTAGE", "Devstage",  "./vst/vst_NR_devstage.png")



### RUVs
differences <- makeGroups(sample_table$Group)
genes <- row.names(raw)
setRUVs <- RUVs(as.matrix(assay(raw)), genes, k=1, differences)
dataRUVs <- DESeqDataSetFromMatrix(countData = setRUVs$normalizedCounts,
colData = sample_table,
design = ~ Group)
RUVs_counts<-t(assay(dataRUVs))
pca_ruvs <- prcomp(RUVs_counts,
center = F,
scale. = F)
pca_data_ruvs <- as.data.frame(pca_ruvs$x)
## RUVs condition
plot_pca(pca_data_ruvs, sample_table$Condition, "RUVs BY CONDITION", pca_ruvs$sdev, "./RUVs/RUVs_by_condition.png", "Condition", sample_table$Group)
# RUVs devstage
plot_pca(pca_data_ruvs, sample_table$DevStage, "RUVs BY DEVSTAGE", pca_ruvs$sdev, "./RUVs/RUVs_by_devstage.png", "Devstage", sample_table$Group)
# RUVs genotype
plot_pca(pca_data_ruvs, sample_table$Genotype, "RUVs BY GENOTYPE", pca_ruvs$sdev, "./RUVs/RUVs_by_genotype.png", "Genotype", sample_table$Group)


## RUVs NR condition
RUVs_NR <- dataRUVs[ ,dataRUVs$Genotype ==  "NR" ]
RUVs_NR_counts<-t(assay(RUVs_NR))
pca_ruvs_NR <- prcomp(RUVs_NR_counts,
center = FALSE,
scale. = FALSE)
pca_data_ruvs_NR <- as.data.frame(pca_ruvs_NR$x)
plot_pca(pca_data_ruvs_NR, RUVs_NR$Condition, "RUVs NR BY CONDITION", pca_ruvs_NR$sdev, "./RUVs/RUVs_NR_condition.png", "Condition", RUVs_NR$Group)
## RUVs NR Devstage
plot_pca(pca_data_ruvs_NR, RUVs_NR$DevStage, "RUVs NR BY DEVSTAGE", pca_ruvs_NR$sdev, "./RUVs/RUVs_NR_devstage.png", "Devstage", RUVs_NR$Group)
## RUVs R devstage
RUVs_R <- dataRUVs[ ,dataRUVs$Genotype ==  "R" ]
RUVs_R_counts<-t(assay(RUVs_R))
pca_ruvs_R <- prcomp(RUVs_R_counts,
center = FALSE,
scale. = FALSE)
pca_data_ruvs_R <- as.data.frame(pca_ruvs_R$x)
plot_pca(pca_data_ruvs_R, RUVs_R$DevStage, "RUVs R BY DEVSTAGE", pca_ruvs_R$sdev, "./RUVs/RUVs_R_devstage.png", "Devstage", RUVs_R$Group)
## RUVs R condition
plot_pca(pca_data_ruvs_R, RUVs_R$Condition, "RUVs R BY CONDITION", pca_ruvs_R$sdev, "./RUVs/RUVs_R_condition.png", "Condition", RUVs_R$Group)


## RUVg
dea <- DESeq(raw)
results_table <- as.data.frame(results(dea))
empirical_control <- rownames(tail(arrange(results_table, pvalue),EC))
RUVgSet <- RUVg(as.matrix(assay(raw)), empirical_control, k = 1)
dataRUVg<- DESeqDataSetFromMatrix(countData = RUVgSet$normalizedCounts,
colData = sample_table,
design = ~ Group)
RUVg_counts<-t(assay(dataRUVg))
pca_ruvg <- prcomp(RUVg_counts,
center = FALSE,
scale. = FALSE)
pca_data_ruvg <- as.data.frame(pca_ruvg$x)
## RUVg devstage
plot_pca(pca_data_ruvg, dataRUVg$DevStage, "RUVg BY DEVSTAGE", pca_ruvg$sdev, "./RUVg/RUVg_devstage.png", "Devstage", dataRUVg$Group)
## RUVg genotype
plot_pca(pca_data_ruvg, dataRUVg$Genotype, "RUVg BY GENOTYPE", pca_ruvg$sdev, "./RUVg/RUVg_genotype.png", "Genotype", dataRUVg$Group)
## RUVg Condition
plot_pca(pca_data_ruvg, dataRUVg$Condition, "RUVg BY CONDITION", pca_ruvg$sdev, "./RUVg/RUVg_condition.png", "Condition", dataRUVg$Group)
## RUVg R condition
RUVg_R <- dataRUVg[ ,dataRUVg$Genotype ==  "R" ]
RUVg_R_counts<-t(assay(RUVg_R))
pca_ruvg_R <- prcomp(RUVg_R_counts,
center = FALSE,
scale. = FALSE)
pca_data_ruvg_R <- as.data.frame(pca_ruvg_R$x)
plot_pca(pca_data_ruvg_R, RUVg_R$Condition, "RUVg R BY CONDITION", pca_ruvg_R$sdev, "./RUVg/RUVg_R_condition.png", "Condition", RUVg_R$Group)
## RUVg R devstage
plot_pca(pca_data_ruvg_R, RUVg_R$DevStage, "RUVg R BY DEVSTAGE", pca_ruvg_R$sdev, "./RUVg/RUVg_R_devstage.png", "Devstage", RUVg_R$Group)
## RUVg NR condition
RUVg_NR <- dataRUVg[ ,dataRUVg$Genotype ==  "NR" ]
RUVg_NR_counts<-t(assay(RUVg_NR))
pca_ruvg_NR <- prcomp(RUVg_NR_counts,
center = FALSE,
scale. = FALSE)
pca_data_ruvg_NR <- as.data.frame(pca_ruvg_NR$x)
plot_pca(pca_data_ruvg_NR, RUVg_NR$Condition, "RUVg NR BY CONDITION", pca_ruvg_NR$sdev, "./RUVg/RUVg_NR_condition.png", "Condition", RUVg_NR$Group)
## RUVg NR devstage
plot_pca(pca_data_ruvg_NR, RUVg_NR$DevStage, "RUVg NR BY DEVSTAGE", pca_ruvg_NR$sdev, "./RUVg/RUVg_NR_devstage.png", "Devstage", RUVg_NR$Group)

## RUVr
design <- model.matrix(~sample_table$Group)
y <- DGEList(counts=counts(raw), group=sample_table$Group)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
res <- residuals(fit, type="deviance")
set_RUVr <- RUVr(y$counts,rownames(y), k=1, res)
dataRUVr <- DESeqDataSetFromMatrix(countData = set_RUVr$normalizedCounts,
colData = sample_table,
design = ~ Group)
RUVr_counts<-t(assay(dataRUVr))
pca_ruvr <- prcomp(RUVr_counts,
center = FALSE,
scale. = FALSE)
pca_data_ruvr <- as.data.frame(pca_ruvr$x)
## RUVr condition
plot_pca(pca_data_ruvr, dataRUVr$Condition, "RUVr BY CONDITION", pca_ruvr$sdev, "./RUVr/RUVr_condition.png", "Condition", dataRUVr$Group)
## RUVr genotype
plot_pca(pca_data_ruvr, dataRUVr$Genotype, "RUVr BY GENOTYPE", pca_ruvr$sdev, "./RUVr/RUVr_genotype.png", "Genotype", dataRUVr$Group)
## RUVr Devstage
plot_pca(pca_data_ruvr, dataRUVr$DevStage, "RUVr BY DEVSTAGE", pca_ruvr$sdev, "./RUVr/RUVr_devstage.png", "Devstagee", dataRUVr$DevStage)

## RUVr R condition
RUVr_R <- dataRUVr[ ,dataRUVr$Genotype ==  "R" ]
RUVr_R_counts<-t(assay(RUVr_R))
pca_ruvr_R <- prcomp(RUVr_R_counts,
center = FALSE,
scale. = FALSE)
pca_data_ruvr_R <- as.data.frame(pca_ruvr_R$x)
plot_pca(pca_data_ruvr_R, RUVr_R$Condition, "RUVr R BY CONDITION", pca_ruvr_R$sdev, "./RUVr/RUVr_R_condition.png", "Condition", RUVr_R$Group)
## RUVr R devstage
plot_pca(pca_data_ruvr_R, RUVr_R$DevStage, "RUVr R BY DEVSTAGE", pca_ruvr_R$sdev, "./RUVr/RUVr_R_devstage.png", "Devstage", RUVr_R$Group)
## RUVr NR condition
RUVr_NR <- dataRUVr[ ,dataRUVr$Genotype ==  "NR" ]
RUVr_NR_counts<-t(assay(RUVr_NR))
pca_ruvr_NR <- prcomp(RUVr_NR_counts,
center = FALSE,
scale. = FALSE)
pca_data_ruvr_NR <- as.data.frame(pca_ruvr_NR$x)
plot_pca(pca_data_ruvr_NR, RUVr_NR$Condition, "RUVr NR BY CONDITION", pca_ruvr_NR$sdev, "./RUVr/RUVr_NR_condition.png", "Condition", RUVr_NR$Group)
## RUVr NR devstage
plot_pca(pca_data_ruvr_NR, RUVr_NR$DevStage, "RUVr NR BY DEVSTAGE", pca_ruvr_NR$sdev, "./RUVr/RUVr_NR_devstage.png", "Devstage", RUVr_NR$Group)

save.image(file='session.Rdata')
