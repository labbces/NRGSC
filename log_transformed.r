# set number of genes for empirical control for removing unwanted variation
EC = 3000
# set working directory
setwd(".")
## install packages
# install.packages(c("BiocManager","readr","dplyr", "magrittr","ggplot2"),
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
library(hexbin, lib.loc =  "./libraries")
library(EDASeq, lib.loc =  "./libraries")
library(RUVSeq, lib.loc =  "./libraries")
library(edgeR, lib.loc = "./libraries")


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

raw_pca <- assay(data)

pca_raw <- prcomp(raw_pca,
                  center = TRUE,
                  scale. = TRUE) 
pca_raw_data <- as.data.frame(pca_raw$rotation)

## plot raw by devstage

p <- ggplot(pca_raw_data, aes(x=PC1, y=PC2, colour=sample_table$Group, shape=sample_table$DevStage)) 
p <- p + geom_point(size=4.5)
p <- p + theme_bw()
p <- p + labs(title = "RAW DATA BY DEVSTAGE", col="Group", shape="DevStage")
p <- p + xlab((paste0("PC1 : ",round(100 * (pca_raw$sdev[1]/sum(pca_raw$sdev)),2),"%")))
p <- p + ylab((paste0("PC1 : ",round(100 * (pca_raw$sdev[2]/sum(pca_raw$sdev)),2),"%")))
p <- p + theme( text = element_text(size=22), 
                panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"))
ggsave( p, filename = "./raw_by_devstage.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)
p
## plot raw Genotype

q <- ggplot(pca_raw_data, aes(x=PC1, y=PC2, colour=sample_table$Group, shape=sample_table$Genotype)) 
q <- q + geom_point(size=4.5)
q <- q + theme_bw()
q <- q + labs(title = "RAW DATA BY Genotype", col="Group", shape="Genotype")
q <- q + xlab((paste0("PC1 : ",round(100 * (pca_raw$sdev[1]/sum(pca_raw$sdev)),2),"%")))
q <- q + ylab((paste0("PC1 : ",round(100 * (pca_raw$sdev[2]/sum(pca_raw$sdev)),2),"%")))
q <- q + theme( text = element_text(size=22), 
                panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"))
q
ggsave( q, filename = "./raw_by_genotype.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

## plot raw condition

r <- ggplot(pca_raw_data, aes(x=PC1, y=PC2, colour=sample_table$Group, shape=sample_table$Condition)) 
r <- r + geom_point(size=4.5)
r <- r + theme_bw()
r <- r + labs(title = "RAW DATA BY CONDITION", col="Group", shape="Condition")
r <- r + xlab((paste0("PC1 : ",round(100 * (pca_raw$sdev[1]/sum(pca_raw$sdev)),2),"%")))
r <- r + ylab((paste0("PC1 : ",round(100 * (pca_raw$sdev[2]/sum(pca_raw$sdev)),2),"%")))
r <- r + theme( text = element_text(size=22), 
                panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"))
r
ggsave( r, filename = "./raw_by_Condition.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)


### pca vst Devstage

pca_vst_devstage <- plotPCA(vst, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes( shape = vst$DevStage))+
  scale_shape_manual(values = seq(0,8)) +
  labs(title = "VST devstage", col="Group", shape="DevStage") +
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
pca_vst_devstage
ggsave( pca_vst_devstage, filename = "./vst/vst_devstage.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

### pca vst genotype

pca_vst_genotype <- plotPCA(vst, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes( shape = vst$Genotype))+
  scale_shape_manual(values = seq(0,8)) +
  labs(title = "VST genotype", col="Group", shape="Genotype") +
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
pca_vst_genotype
ggsave( pca_vst_genotype, filename = "./vst/vst_genotype.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

### pca vst condition

pca_vst_condition <- plotPCA(vst, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes( shape = vst$Condition))+
  scale_shape_manual(values = seq(0,8)) +
  labs(title = "VST genotype", col="Group", shape="Condition") +
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
pca_vst_condition
ggsave( pca_vst_genotype, filename = "./vst/vst_condition.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

### sub plot vst R condition
vst_R_condition <- vst[ ,vst$Genotype ==  "R" ]
pca_vst_R_condition <- plotPCA(vst_R_condition, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes( shape = vst_R_condition$Condition))+
  scale_shape_manual(values = seq(0,8)) +
  labs(title = "VST R CONDITION", col="Group", shape="Condition") +
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
pca_vst_R_condition
ggsave(pca_vst_R_condition, filename = "./vst/vst_R_condition.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

### sub plot vst R devstage
VST_R <- vst[ ,vst$Genotype ==  "R" ]
pca_vst_R_devstage <- plotPCA(VST_R, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes( shape = VST_R$DevStage))+
  scale_shape_manual(values = seq(0,8)) +
  labs(title = "VST R Devstage", col="Group", shape="Devstage") +
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
pca_vst_R_devstage
ggsave( pca_vst_R_devstage, filename = "./vst/vst_R_devstage.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)


### subplot vst NR devstage

VST_NR <- vst[ ,vst$Genotype ==  "NR" ]
pca_vst_NR_devstage <- plotPCA(VST_NR, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes( shape = VST_NR$DevStage))+
  scale_shape_manual(values = seq(0,8)) +
  labs(title = "VST NR DEVSTAGE", col="Group", shape="DevStage") +
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
pca_vst_NR_devstage
ggsave( pca_vst_NR_devstage, filename = "./vst/vst_devstage_NR.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

### subplot vst NR condition

pca_vst_NR_condition <- plotPCA(VST_NR, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes( shape = VST_NR$Condition))+
  scale_shape_manual(values = seq(0,8)) +
  labs(title = "VST NR CONDITION", col="Group", shape="Condition") +
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
pca_vst_NR_condition
ggsave( pca_vst_NR_condition, filename = "./vst/vst_NR_condition.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)


differences <- makeGroups(sample_table$Group)
genes <- row.names(raw)
setRUVs <- RUVs(as.matrix(assay(raw)), genes, k=1, differences)
dataRUVs <- DESeqDataSetFromMatrix(countData = setRUVs$normalizedCounts,
                                   colData = sample_table,
                                   design = ~ Group)

dataRUVs <- estimateSizeFactors(dataRUVs)
ruvs_vst <- varianceStabilizingTransformation(dataRUVs)

## RUVs condition

pca_ruvs_condition <- plotPCA(ruvs_vst, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes(shape = ruvs_vst$Condition, color = ruvs_vst$Group))+
  scale_shape_manual(values=seq(0,8)) +
  labs(title = "RUVs VST condition", col="Group", shape="Condition")+
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )

pca_ruvs_condition
ggsave( pca_ruvs_condition, filename = "./RUVs/RUVs_vst_condition.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

## RUVs genotype
pca_ruvs_genotype <- plotPCA(ruvs_vst, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes(shape = ruvs_vst$Genotype, color = ruvs_vst$Group))+
  scale_shape_manual(values=seq(0,8)) +
  labs(title = "RUVs VST genotype", col="Group", shape="Genotype")+
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )

pca_ruvs_genotype
ggsave( pca_ruvs_genotype, filename = "./RUVs/RUVs_vst_genotype.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

## RUVs Devstage

pca_ruvs_devstage <- plotPCA(ruvs_vst, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes(shape = ruvs_vst$DevStage, color = ruvs_vst$Group))+
  scale_shape_manual(values=seq(0,8)) +
  labs(title = "RUVs VST Devstage", col="Group", shape="Devstage")+
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )

pca_ruvs_devstage
ggsave( pca_ruvs_devstage, filename = "./RUVs/RUVs_vst_devstage.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

## RUVs NR condition

RUVs_NR <- ruvs_vst[ ,ruvs_vst$Genotype ==  "NR" ]
ruvs_NR_condition <- plotPCA(RUVs_NR, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes( shape = RUVs_NR$Condition))+
  scale_shape_manual(values = seq(0,8)) +
  labs(title = "RUVs NR CONDITION", col="Group", shape="Condition") +
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
ruvs_NR_condition
ggsave( ruvs_NR_condition, filename = "./RUVs/NR/RUVs_NR_condition.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

## RUVs R devstage
RUVs_R <- ruvs_vst[ ,ruvs_vst$Genotype ==  "R" ]
ruvs_R_devstage <- plotPCA(RUVs_R, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes( shape = RUVs_R$DevStage))+
  scale_shape_manual(values = seq(0,8)) +
  labs(title = "RUVs R Devstage", col="Group", shape="Devstage") +
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
ruvs_R_devstage
ggsave( ruvs_R_devstage, filename = "./RUVs/R/RUVs_R_devstage.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

## RUVs R condition
RUVs_R <- ruvs_vst[ ,ruvs_vst$Genotype ==  "R" ]
ruvs_R_condition <- plotPCA(RUVs_R, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes( shape = RUVs_R$Condition))+
  scale_shape_manual(values = seq(0,8)) +
  labs(title = "RUVs R CONDITION", col="Group", shape="Condition") +
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
ruvs_R_condition
ggsave( ruvs_R_condition, filename = "./RUVs/R/RUVs_R_condition.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

## RUVs NR devstage

RUVs_NR <- ruvs_vst[ ,ruvs_vst$Genotype ==  "NR" ]
ruvs_NR_devstage <- plotPCA(RUVs_R, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes( shape = RUVs_NR$DevStage))+
  scale_shape_manual(values = seq(0,8)) +
  labs(title = "RUVs NR Devstage", col="Group", shape="Devstage") +
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
ruvs_NR_devstage
ggsave( ruvs_NR_devstage, filename = "./RUVs/NR/RUVs_NR_devstage.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

dea <- DESeq(raw)
results_table <- as.data.frame(results(dea))
empirical_control <- rownames(tail(arrange(results_table, pvalue),EC))
RUVgSet <- RUVg(as.matrix(assay(raw)), empirical_control, k = 1)
dataRUVg<- DESeqDataSetFromMatrix(countData = RUVgSet$normalizedCounts,
                                  colData = sample_table,
                                  design = ~ Group)

dataRUVg <- estimateSizeFactors(dataRUVg)
ruvg_vst <- varianceStabilizingTransformation(dataRUVg)

## RUVg devstage

ruvg_devstage <- plotPCA(ruvg_vst, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes(shape = ruvg_vst$DevStage, color = ruvg_vst$Group))+
  scale_shape_manual(values=seq(0,8)) +
  labs(title = "RUVg Devstage", col="Group", shape="DevStage")+
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )

ruvg_devstage
ggsave( ruvg_devstage, filename = "./RUVg/ruvg_devstage.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)


## RUVg genotype

ruvg_genotype <- plotPCA(ruvg_vst, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes(shape = ruvg_vst$Genotype, color = ruvg_vst$Group))+
  scale_shape_manual(values=seq(0,8)) +
  labs(title = "RUVg genotype", col="Group", shape="Genotype")+
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )

ruvg_genotype
ggsave( ruvg_genotype, filename = "./RUVg/ruvg_genotype.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)


## RUVg Condition

ruvg_condition <- plotPCA(ruvg_vst, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes(shape = ruvg_vst$Condition, color = ruvg_vst$Group))+
  scale_shape_manual(values=seq(0,8)) +
  labs(title = "RUVg Condition", col="Group", shape="Condition")+
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )

ruvg_condition
ggsave( ruvg_condition, filename = "./RUVg/ruvg_condition.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

## RUVg R condition

RUVg_R <- ruvg_vst[ ,ruvg_vst$Genotype ==  "R" ]
ruvg_R_condition <- plotPCA(RUVg_R, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes( shape = RUVg_R$Condition))+
  scale_shape_manual(values = seq(0,8)) +
  labs(title = "RUVg R condition", col="Group", shape="Condition") +
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
ruvg_R_condition
ggsave( ruvg_R_condition, filename = "./RUVg/R/ruvg_R_condition.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)


## RUVg R devstage

RUVg_R <- ruvg_vst[ ,ruvg_vst$Genotype ==  "R" ]
ruvg_R_devstage <- plotPCA(RUVg_R, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes( shape = RUVg_R$DevStage))+
  scale_shape_manual(values = seq(0,8)) +
  labs(title = "RUVg R devstage", col="Group", shape="DevStage") +
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
ruvg_R_devstage
ggsave( ruvg_R_devstage, filename = "./RUVg/R/ruvg_R_devstage.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

## RUVg NR condition 
RUVg_NR <- ruvg_vst[ ,ruvg_vst$Genotype ==  "NR" ]
ruvg_NR_condition <- plotPCA(RUVg_NR, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes( shape = RUVg_NR$Condition))+
  scale_shape_manual(values = seq(0,8)) +
  labs(title = "RUVg NR condition", col="Group", shape="Condition") +
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
ruvg_NR_condition
ggsave( ruvg_NR_condition, filename = "./RUVg/NR/ruvg_NR_condition.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

## RUVg NR devstage

RUVg_NR <- ruvg_vst[ ,ruvg_vst$Genotype ==  "NR" ]
ruvg_NR_devstage <- plotPCA(RUVg_NR, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes( shape = RUVg_NR$DevStage))+
  scale_shape_manual(values = seq(0,8)) +
  labs(title = "RUVg NR devstage", col="Group", shape="DevStage") +
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
ruvg_NR_devstage
ggsave( ruvg_NR_devstage, filename = "./RUVg/NR/ruvg_NR_devstage.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)
##
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
RUVr_vst <- varianceStabilizingTransformation(dataRUVr)

## RUVr condition
ruvr_condition <- plotPCA(RUVr_vst, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes(shape = RUVr_vst$Condition))+
  scale_shape_manual(values=seq(0,8)) +
  labs(title = "RUVr Condition", col="Group", shape="Condition")+
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
ruvr_condition
ggsave(ruvr_condition, filename = "./RUVr/RUVr_condition.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)
## RUVr genotype
ruvr_genotype <- plotPCA(RUVr_vst, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes(shape = RUVr_vst$Genotype))+
  scale_shape_manual(values=seq(0,8)) +
  labs(title = "RUVr Genotype", col="Group", shape="Genotype")+
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
ruvr_genotype
ggsave(ruvr_condition, filename = "./RUVr/RUVr_genotype.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

## RUVr Devstage
ruvr_devstage <- plotPCA(RUVr_vst, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes(shape = RUVr_vst$DevStage))+
  scale_shape_manual(values=seq(0,8)) +
  labs(title = "RUVr Devstage", col="Group", shape="Devstage")+
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
ruvr_devstage
ggsave(ruvr_condition, filename = "./RUVr/RUVr_devstage.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

## RUVr R condition
RUVr_R <- RUVr_vst[ ,RUVr_vst$Genotype ==  "R" ]
ruvr_R_condition <- plotPCA(RUVr_R, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes(shape = RUVr_R$Condition))+
  scale_shape_manual(values=seq(0,8)) +
  labs(title = "RUVr R CONDITION", col="Group", shape="Condition")+
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
ruvr_R_condition
ggsave(ruvr_R_condition, filename = "./RUVr/R/RUVr_R_condition.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)


## RUVr R devstage
RUVr_R <- RUVr_vst[ ,RUVr_vst$Genotype ==  "R" ]
ruvr_R_devstage <- plotPCA(RUVr_R, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes(shape = RUVr_R$DevStage))+
  scale_shape_manual(values=seq(0,8)) +
  labs(title = "RUVr R DEVSTAGE", col="Group", shape="Devstage")+
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
ruvr_R_devstage
ggsave(ruvr_R_devstage, filename = "./RUVr/R/RUVr_R_devstage.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

## RUVr NR condition 
RUVr_NR <- RUVr_vst[ ,RUVr_vst$Genotype ==  "NR" ]
ruvr_NR_condition <- plotPCA(RUVr_NR, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes(shape = RUVr_NR$Condition))+
  scale_shape_manual(values=seq(0,8)) +
  labs(title = "RUVr NR CONDITION", col="Group", shape="Condition")+
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
ruvr_NR_condition
ggsave(ruvr_NR_condition, filename = "./RUVr/NR/RUVr_NR_condition.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

## RUVr NR devstage

RUVr_NR <- RUVr_vst[ ,RUVr_vst$Genotype ==  "NR" ]
ruvr_NR_devstage <- plotPCA(RUVr_NR, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes(shape = RUVr_NR$DevStage))+
  scale_shape_manual(values=seq(0,8)) +
  labs(title = "RUVr NR DEVSTAGE", col="Group", shape="Devstage")+
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
ruvr_NR_devstage
ggsave(ruvr_NR_devstage, filename = "./RUVr/NR/RUVr_NR_devstage.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

