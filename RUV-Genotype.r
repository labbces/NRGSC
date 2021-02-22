<<<<<<< HEAD
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
  BiocManager::install("edgeR", lib = "./libraries")

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
rld <- rlog(data, blind = FALSE)

### pca vst
pca_vst <- plotPCA(vst, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes( shape = vst$DevStage))+
  scale_shape_manual(values = seq(0,8)) +
  labs(title = "VST", col="Group", shape="DevStage") +
    theme( text = element_text(size=22), 
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
 	  panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.line = element_line(colour = "black") )
pca_vst
ggsave( pca_vst, filename = "./vst/vst.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

### sub plot genotype R
VST_R <- vst[ ,vst$Genotype ==  "R" ]
pca_vst_R <- plotPCA(VST_R, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes( shape = VST_R$DevStage))+
  scale_shape_manual(values = seq(0,8)) +
  labs(title = "VST GENOTYPE R", col="Group", shape="DevStage") +
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
pca_vst_R
ggsave( pca_vst_R, filename = "./vst/vst_R.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

### subplot genotype vst NR

VST_NR <- vst[ ,vst$Genotype ==  "NR" ]
pca_vst_NR <- plotPCA(VST_NR, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes( shape = VST_NR$DevStage))+
  scale_shape_manual(values = seq(0,8)) +
  labs(title = "VST GENOTYPE NR", col="Group", shape="DevStage") +
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
pca_vst_NR
ggsave( pca_vst_R, filename = "./vst/vst_NR.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)


#### pcs rld
pca_rld <- plotPCA(rld, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes(shape = rld$DevStage)) +
  scale_shape_manual(values=seq(0,8)) +
  labs(title = "RLD", col = "Group", shape = "DevStage")+
  theme_bw() +
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
pca_rld
ggsave( pca_rld, filename = "./rld/rld.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

#### subplots rld genotype N
RLD_R <- rld[ ,rld$Genotype ==  "R" ]
pca_rld_R <- plotPCA(RLD_R, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes(shape = RLD_R$DevStage)) +
  scale_shape_manual(values=seq(0,8)) +
  labs(title = "RLD GENOTYPE R", col = "Group", shape = "DevStage")+
  theme_bw() +
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
pca_rld_R
ggsave( pca_rld_R, filename = "./rld/rld_R.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

### subplot rld genotype NR

RLD_NR <- rld[ ,rld$Genotype ==  "NR" ]
pca_rld_NR <- plotPCA(RLD_NR, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes(shape = RLD_NR$DevStage)) +
  scale_shape_manual(values=seq(0,8)) +
  labs(title = "RLD GENOTYPE NR", col = "Group", shape = "DevStage")+
  theme_bw() +
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
pca_rld_NR
ggsave( pca_rld_NR, filename = "./rld/rld_NR.png",units = "cm",width = 15*1.3, height = 15, dpi = 320)


#### RUVg
dea <- DESeq(raw)
results_table <- as.data.frame(results(dea))
empirical_control <- rownames(tail(arrange(results_table, pvalue),EC))
RUVgSet <- RUVg(as.matrix(assay(raw)), empirical_control, k = 1)
dataRUVg<- DESeqDataSetFromMatrix(countData = RUVgSet$normalizedCounts,
                                 colData = sample_table,
                                 design = ~ Group)

dataRUVg <- estimateSizeFactors(dataRUVg)
ruvg_vst <- varianceStabilizingTransformation(dataRUVg)
ruvg_rld <- rlog(dataRUVg, blind = FALSE)
  
# PCA RUVg VST
pca_ruv_vst <- plotPCA(ruvg_vst, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes(shape = ruvg_vst$DevStage, color = ruvg_vst$Group))+
  scale_shape_manual(values=seq(0,8)) +
  labs(title = "RUVg VST", col="Group", shape="DevStage")+
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )

pca_ruv_vst
ggsave( pca_ruv_vst, filename = "./RUVg/RUVg_vst.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

## subplot RUVg VST genotype NR

VST_RUVg_NR <- ruvg_vst[ ,ruvg_vst$Genotype ==  "NR" ]
pca_ruvg_vst_NR <- plotPCA(VST_RUVg_NR, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes( shape = VST_RUVg_NR$DevStage))+
  scale_shape_manual(values = seq(0,8)) +
  labs(title = "RUVg VST GENOTYPE NR", col="Group", shape="DevStage") +
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
pca_ruvg_vst_NR
ggsave( pca_ruvg_vst_NR, filename = "./RUVg/NR/RUVg_vst_NR.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

## subplot RUVg VST genotype R

VST_RUVg_R <- ruvg_vst[ ,ruvg_vst$Genotype ==  "R" ]
pca_ruvg_vst_R <- plotPCA(VST_RUVg_R, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes( shape = VST_RUVg_R$DevStage))+
  scale_shape_manual(values = seq(0,8)) +
  labs(title = "RUVg VST GENOTYPE R", col="Group", shape="DevStage") +
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
pca_ruvg_vst_R
ggsave( pca_ruvg_vst_NR, filename = "./RUVg/R/RUVg_vst_R.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

# PCA RUVg RLD
pca_ruvg_rld <- plotPCA(ruvg_rld, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes(shape = ruvg_rld$DevStage))+
  scale_shape_manual(values=seq(0,8)) +
  labs(title = "RUVg RLD", col="Group", shape="DevStage")+
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
pca_ruvg_rld
ggsave(pca_ruvg_rld, filename = "./RUVg/RUVg_rld.png",units = "cm",width = 15*1.3, height = 15, dpi = 320)

## RUVs

differences <- makeGroups(sample_table$Group)
genes <- row.names(raw)
setRUVs <- RUVs(as.matrix(assay(raw)), genes, k=1, differences)
dataRUVs <- DESeqDataSetFromMatrix(countData = setRUVs$normalizedCounts,
                                  colData = sample_table,
                                  design = ~ Group)

dataRUVs <- estimateSizeFactors(dataRUVs)
ruvs_vst <- varianceStabilizingTransformation(dataRUVs)
ruvs_rld <- rlog(dataRUVs, blind = FALSE)

# PCA RUVs VST
pca_ruvs_vst <- plotPCA(ruvs_vst, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes(shape = ruvs_vst$DevStage, color = ruvs_vst$Group))+
  scale_shape_manual(values=seq(0,8)) +
  labs(title = "RUVs VST", col="Group", shape="DevStage")+
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )

pca_ruvs_vst
ggsave( pca_ruvs_vst, filename = "./RUVs/RUVs-vst.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

## subplot PCA RUVs VST genotype NR

VST_RUVs_NR <- ruvs_vst[ ,ruvs_vst$Genotype ==  "NR" ]
pca_ruvs_vst_NR <- plotPCA(VST_RUVs_NR, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes( shape = VST_RUVs_NR$DevStage))+
  scale_shape_manual(values = seq(0,8)) +
  labs(title = "RUVs VST GENOTYPE NR", col="Group", shape="DevStage") +
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
pca_ruvs_vst_NR
ggsave( pca_ruvs_vst_NR, filename = "./RUVs/NR/RUVs_vst_NR.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

## subplot PCA RUVs VST genotype R

VST_RUVs_R <- ruvs_vst[ ,ruvs_vst$Genotype ==  "R" ]
pca_ruvs_vst_R <- plotPCA(VST_RUVs_R, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes( shape = VST_RUVs_R$DevStage))+
  scale_shape_manual(values = seq(0,8)) +
  labs(title = "RUVs VST GENOTYPE R", col="Group", shape="DevStage") +
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
pca_ruvs_vst_R
ggsave( pca_ruvs_vst_R, filename = "./RUVs/R/RUVs_vst_R.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)


# PCA RUVs RLD
pca_ruvs_rld <- plotPCA(ruvs_rld, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes(shape = ruvs_rld$DevStage))+
  scale_shape_manual(values=seq(0,8)) +
  labs(title = "RUVs RLD", col="Group", shape="DevStage")+
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
pca_ruvs_rld
ggsave(pca_ruvs_rld, filename = "./RUVs/RUVs_rld.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

# RUVr residuals

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
RUVr_rld <- rlog(dataRUVr, blind = FALSE)

## graph vst RUVr
pca_ruvr_vst <- plotPCA(RUVr_vst, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes(shape = RUVr_vst$DevStage))+
  scale_shape_manual(values=seq(0,8)) +
  labs(title = "RUVr VST", col="Group", shape="DevStage")+
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
pca_ruvr_vst
ggsave(pca_ruvr_vst, filename = "./RUVr/RUVr_vst.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

## graph rld RUVr

pca_ruvr_rld <- plotPCA(RUVr_rld, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes(shape = RUVr_rld$DevStage))+
  scale_shape_manual(values=seq(0,8)) +
  labs(title = "RUVr RLD", col="Group", shape="DevStage")+
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
pca_ruvr_rld
ggsave(pca_ruvr_rld, filename = "./RUVr/RUVr_rld.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

## graph rld RUVr genotype R

RUVr_rld_R <- RUVr_rld[ ,RUVr_rld$Genotype ==  "R" ]
pca_ruvr_rld_R <- plotPCA(RUVr_rld_R, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes( shape = RUVr_rld_R$DevStage))+
  scale_shape_manual(values = seq(0,8)) +
  labs(title = "RUVr RLD GENOTYPE R", col="Group", shape="DevStage") +
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
pca_ruvr_rld_R
ggsave( pca_ruvr_rld_R, filename = "./RUVr/R/RUVr_rld_R.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

## graph rld RUVr genotype NR

RUVr_rld_NR <- RUVr_rld[ ,RUVr_rld$Genotype ==  "NR" ]
pca_ruvr_rld_NR <- plotPCA(RUVr_rld_NR, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes( shape = RUVr_rld_NR$DevStage))+
  scale_shape_manual(values = seq(0,8)) +
  labs(title = "RUVr RLD GENOTYPE NR", col="Group", shape="DevStage") +
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
pca_ruvr_rld_NR
ggsave( pca_ruvr_rld_R, filename = "./RUVr/NR/RUVr_rld_NR.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

## graph vst RUVr genotype NR
RUVr_vst_NR <- RUVr_vst[ ,RUVr_vst$Genotype ==  "NR" ]
pca_ruvr_vst_NR <- plotPCA(RUVr_vst_NR, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes(shape = RUVr_vst_NR$DevStage))+
  scale_shape_manual(values=seq(0,8)) +
  labs(title = "RUVr VST", col="Group", shape="DevStage")+
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
pca_ruvr_vst_NR
ggsave(pca_ruvr_vst_NR, filename = "./RUVr/NR/RUVr_vst_NR.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

## graph vst RUVr genotype R
RUVr_vst_R <- RUVr_vst[ ,RUVr_vst$Genotype ==  "R" ]
pca_ruvr_vst_R <- plotPCA(RUVr_vst_R, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes(shape = RUVr_vst_R$DevStage))+
  scale_shape_manual(values=seq(0,8)) +
  labs(title = "RUVr VST", col="Group", shape="DevStage")+
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
pca_ruvr_vst_R
ggsave(pca_ruvr_vst_R, filename = "./RUVr/R/RUVr_vst_R.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

=======
# set working directory
setwd("C://Users/Avell 5/Downloads/NRGSC-main/NRGSC-main/")
## install packages
# install.packages(c("BiocManager","readr","dplyr", "magrittr","ggplot2","hexbin"),
# lib = "./libraries",
# repos = 'http://cran.us.r-project.org')
 
# BiocManager::install("tximport", lib = "./libraries")
# BiocManager::install("DESeq2", lib = "./libraries")
# BiocManager::install("RUVSeq", lib = "./libraries")
# BiocManager::install("EDASeq", lib = "./libraries")

#  load libraries
# library(readr, lib.loc =  "./libraries")
# library(dplyr, lib.loc =  "./libraries/")
# library(magrittr,  lib.loc = "./libraries")
# library(tximport, lib.loc = "./libraries")
# library(DESeq2, lib.loc = "./libraries")
# library(ggplot2, lib.loc =  "./libraries")
# library(hexbin, lib.loc =  "./libraries")
# library(EDASeq, lib.loc =  "./libraries")
# library(RUVSeq, lib.loc =  "./libraries")
 
 library(readr)
 library(dplyr)
 library(magrittr)
 library(tximport)
 library(DESeq2)
 library(ggplot2)
 library(hexbin)
 library(EDASeq)
 library(RUVSeq)

# load metadata
sample_table <-read.table("metadata_complete.csv", sep = ",", header = T)
# delete_samples
#sample_table<-rows_delete(sample_table, tibble(Sample.Number = 37))
# load files paths
sample_files = paste0(pull(sample_table , "Sample_file"), "/quant.sf")
# name table columns
names(sample_files) = pull(sample_table, "Sample_file")
# relate genes to transcripts
tx2gene = read.table("txgen2.txt", sep = "\t", col.names =c("genid","transid"))
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
rld <- rlog(data, blind = FALSE)

### pca vst
pca_vst <- plotPCA(vst, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes(shape=vst$DevStage))
#  geom_point(size=4.5, aes(colour = vst$DevStage, 
#                           shape = vst$Genotype))+
#    theme( text = element_text(size=22), 
#    panel.border = element_blank(),
#    panel.grid.major = element_blank(),
# 	  panel.grid.minor = element_blank(),
#    plot.title = element_text(hjust = 0.5),
#    axis.line = element_line(colour = "black") )

ggsave( pca_vst, filename = "./PCA_VST.png",units = "cm",width = 20*1.3, height = 20,dpi = 320)

#### pcs rld

pca_rld <- plotPCA(rld, intgroup ='Group') +
  theme_bw() +
  geom_point(size=4.5, aes(shape=rld$DevStage))
#pca_rld <- plotPCA(rld, intgroup ='Group') +
#  geom_point(size=4.5, aes(colour = rld$Condition, shape = rld$Genotype))+
#  labs(title = "RLD", col="Condition", shape="Genotype")+
#  theme_bw() +
#  theme( text = element_text(size=22), 
#         panel.border = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.line = element_line(colour = "black") )

ggsave( pca_rld, filename = "./PCA_nrld.png",units = "cm",width = 20*1.3, height = 20,dpi = 320)


#### RUV

differences <- makeGroups(sample_table$Condition)
genes <- row.names(count_data$counts)
setRUV <- RUVs(as.matrix(assay(raw)), genes, k=1, differences)
dataRUV <- DESeqDataSetFromMatrix(countData = setRUV$normalizedCounts,
                                 colData = sample_table,
                                 design = ~ Condition)

dataRUV <- estimateSizeFactors(dataRUV)
ruv_vst <- varianceStabilizingTransformation(dataRUV)
ruv_rld <- rlog(dataRUV, blind = FALSE)
  

# PCA RUV VST
pca_ruv_vst <- plotPCA(ruv_vst, intgroup ='Condition') +
  theme_bw() +
  geom_point(size=4.5, aes(colour = Condition, shape = ruv_vst$Genotype))+
  labs(title = "PCA RUV CONDITION VST", col="Condition", shape="Genotype")+
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )

pca_ruv_vst
ggsave( pca_ruv_vst, filename = "./PCA_RUV-CONDITION-VST.png",units = "cm",width = 20*1.3, height = 20,dpi = 320)

# PCA RUV RLD

pca_ruv_rld <- plotPCA(ruv_rld, intgroup ='Condition') +
  theme_bw() +
  geom_point(size=4.5, aes(colour = Condition, shape = ruv_rld$Genotype))+
  labs(title = "RUV CONDITION RLD", col="Condition", shape="Genotype")+
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )

ggsave( pca_ruv_rld, filename = "./PCA_RUV-CONDITION_RLD.png",units = "cm",width = 20*1.3, height = 20,dpi = 320)

## RUV GENOTYPE

differences <- makeGroups(sample_table$Genotype)
genes <- row.names(count_data$counts)
setRUV <- RUVs(as.matrix(assay(raw)), genes, k=1, differences)
dataRUV <- DESeqDataSetFromMatrix(countData = setRUV$normalizedCounts,
                                  colData = sample_table,
                                  design = ~ Condition)
dataRUV <- estimateSizeFactors(dataRUV)
ruv_vst <- varianceStabilizingTransformation(dataRUV)
ruv_rld <- rlog(dataRUV, blind = FALSE)

# PCA RUV RLD GENOTYPE
pca_ruv_rld <- plotPCA(ruv_rld, intgroup ='Condition') +
  theme_bw() +
  geom_point(size=4.5, aes(colour = Condition, shape = ruv_rld$Genotype))+
  labs(title = "ruv rld", col="Condition", shape="Genotype")+
  theme( text = element_text(size=22),
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )
ggsave( pca_ruv_rld, filename = "./PCA_RUV_RLD-GENOTYPE.png",units = "cm",width = 20*1.3, height = 20,dpi = 320)

# PCA RUV VST GENOTYPE
pca_ruv_vst <- plotPCA(ruv_vst, intgroup ='Condition') +
  theme_bw() +
  geom_point(size=4.5, aes(colour = Condition, shape = ruv_vst$Genotype))+
  labs(title = "PCA RUV~CONDITION VST", col="Condition", shape="Genotype")+
  theme( text = element_text(size=22),
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )

ggsave( pca_ruv_vst, filename = "./PCA_RUV-VST-GENOTYPE.png",units = "cm",width = 20*1.3, height = 20,dpi = 320)
>>>>>>> 1c1b9ae06c0dca3c986ea0ab107d4e153bae9594
