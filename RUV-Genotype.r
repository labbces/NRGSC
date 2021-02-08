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
pca_vst <- plotPCA(vst, intgroup ='Condition') +
  theme_bw() +
  geom_point(size=4.5, aes(colour = Condition, shape = vst$Genotype))+
  labs(title = "VST", col="Condition", shape="Genotype")+
  theme( text = element_text(size=22), 
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
 	  panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.line = element_line(colour = "black") )

ggsave( pca_vst, filename = "./PCA_VST.png",units = "cm",width = 20*1.3, height = 20,dpi = 320)

#### pcs rld

pca_rld <- plotPCA(rld, intgroup ='Condition') +
  geom_point(size=4.5, aes(colour = Condition, shape = rld$Genotype))+
  labs(title = "RLD", col="Condition", shape="Genotype")+
  theme_bw() +
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.line = element_line(colour = "black") )

pca_rld
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
