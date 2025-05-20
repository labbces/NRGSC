## variables
wd = "/Storage/data1/jorge.munoz/NRGSC.new/code"
cores = 12
# set working directory
setwd(wd)

## install packages
# install.packages(c("BiocManager","backports","tidyverse"),
# lib = "./libraries",
# repos = 'http://cran.us.r-project.org')
#BiocManager::install("tximport", lib = "./../libraries")
#BiocManager::install("DESeq2", lib = "./../libraries")
#BiocManager::install("BiocParallel", lib = "./../libraries")

## load libraries 
library(tximport, lib.loc = "./../libraries")
library(DESeq2, lib.loc = "./../libraries")
library(backports, lib.loc = "./../libraries")
library(tidyverse, lib.loc = "./../libraries")
library(BiocParallel, lib.loc = "./../libraries")

# parallel envirorment
register(MulticoreParam(cores))

# define fuction to make DEA analysis in DESeq2
get_full_gene_count_matrixes <- function(metadata, txgen2) {
        sample_table <-read.table(metadata, sep = ",", header = T)
        sample_table$Group <- as.factor(paste(sample_table$Genotype, '_', sample_table$Condition, '_', sample_table$DevStage, sep='')) 
        # load files paths
        sample_files = paste0("/Storage/data1/jorge.munoz/NRGSC/", pull(sample_table , "Sample_file"), "/quant.sf")
        # name table columns
        names(sample_files) = pull(sample_table, "Sample_file")
        # relate genes to transcripts
        tx2gene = read.table(txgen2, sep = ",", col.names =c("transid","geneid"))
        # import count data to tximport
        count_data = tximport( files = sample_files,
                               type = "salmon",
                               tx2gene =  tx2gene,
			       ignoreTxVersion = F)	

        sample_table$Condition <- as.factor((sample_table$Condition))
        sample_table$Genotype <- as.factor((sample_table$Genotype))
        sample_table$DevStage <- as.factor((sample_table$DevStage))
        sample_table$Individual <- as.factor((sample_table$Individual))
        raw <- DESeqDataSetFromTximport(txi = count_data,
                                        colData = sample_table,
                                        design = ~ Group)
        # full raw and normalized count matrixes for all experiments
        data <- estimateSizeFactors(raw)
        vst <- varianceStabilizingTransformation(data)                	
	df_data <- assay(raw)
	df_vst <- assay(vst)
	write.table(df_data, file = paste0("./../results/gene/", "all_raw_counts.tsv", sep = ""), sep = "\t" ,quote = F, col.names = T)
        write.table(df_vst, file = paste0("./../results/gene/", "all_vst_counts.tsv", sep = ""), sep = "\t" ,quote = F, col.names = T)       
}

get_full_trans_count_matrixes <- function(metadata, txgen2) {
        sample_table <-read.table(metadata, sep = ",", header = T)
        sample_table$Group <- as.factor(paste(sample_table$Genotype, '_', sample_table$Condition, '_', sample_table$DevStage, sep=''))
        # load files paths
        sample_files = paste0("/Storage/data1/jorge.munoz/NRGSC/", pull(sample_table , "Sample_file"), "/quant.sf")
        # name table columns
        names(sample_files) = pull(sample_table, "Sample_file")
        # relate genes to transcripts
        tx2gene = read.table(txgen2, sep = ",", col.names =c("transid","geneid"))
        # import count data to tximport
        count_data = tximport( files = sample_files,
                               type = "salmon",
                               tx2gene =  tx2gene,
			       txOut = T,	
                               ignoreTxVersion = F)
        sample_table$Condition <- as.factor((sample_table$Condition))
        sample_table$Genotype <- as.factor((sample_table$Genotype))
        sample_table$DevStage <- as.factor((sample_table$DevStage))
        sample_table$Individual <- as.factor((sample_table$Individual))
        raw <- DESeqDataSetFromTximport(txi = count_data,
                                        colData = sample_table,
                                        design = ~ Group)
        # full raw and normalized count matrixes for all experiments
        data <- estimateSizeFactors(raw)
        vst <- varianceStabilizingTransformation(data)
        df_data <- assay(raw)
        df_vst <- assay(vst)
        write.table(df_data, file = paste0("./../results/trans/", "all_raw_counts.tsv", sep = ""), sep = "\t" ,quote = F, col.names = T)
        write.table(df_vst, file = paste0("./../results/trans/", "all_vst_counts.tsv", sep = ""), sep = "\t" ,quote = F, col.names = T)
}
get_full_gene_count_matrixes("./../data/metadata_complete.csv",  "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")
get_full_trans_count_matrixes("./../data/metadata_complete.csv",  "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")
