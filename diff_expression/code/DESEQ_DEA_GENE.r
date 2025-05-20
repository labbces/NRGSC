## variables
wd = "/Storage/data1/jorge.munoz/NRGSC.new/code"
outdir = "/Storage/data1/jorge.munoz/NRGSC.new/perGENOTYPE/gene/"
cores = 6
# set working directory
setwd(wd)
# make results dir
dir.create(outdir)
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
dea_by_group <- function(metadata, group1, group2, txgen2) {
	dir.create(paste0(outdir, group1, "_vs_", group2))
        sample_table <-read.table(metadata, sep = ",", header = T)
        sample_table$Group <- as.factor(paste(sample_table$Genotype, '_', sample_table$Condition, '_', sample_table$DevStage, sep=''))
        
        filter <- sample_table %>% filter( Group == group1 | Group == group2 )

        # load files paths
        sample_files = paste0("/Storage/data1/jorge.munoz/NRGSC/", pull(filter , "Sample_file"), "/quant.sf")
        # name table columns
        names(sample_files) = pull(filter, "Sample_file")
        # relate genes to transcripts
        tx2gene = read.table(txgen2, sep = ",", col.names =c("transid","geneid"))
	# GENE MODE 
        # import count data to tximport
        count_data = tximport( files = sample_files,
                               type = "salmon",
                               tx2gene =  tx2gene,
                               ignoreTxVersion = F,
			       ignoreAfterBar = T)
        
        sample_table$Condition <- as.factor((sample_table$Condition))
        sample_table$Genotype <- as.factor((sample_table$Genotype))
        sample_table$DevStage <- as.factor((sample_table$DevStage))
        sample_table$Individual <- as.factor((sample_table$Individual))
        
        raw <- DESeqDataSetFromTximport(txi = count_data,
                                        colData = filter,
                                        design = ~ Group)
        dim(raw)
        temp <- as.data.frame(counts(raw))
        logic <-(apply(temp,c(1,2), function(x){x>0}))
        filter_genes <- rowSums(logic)>3
        fi <- raw[filter_genes,]
        dim(fi)
        temp <- NULL
        data <- estimateSizeFactors(fi)
        vst <- varianceStabilizingTransformation(data)
        
	# regex names
	# g1-g2-filtered_raw_counts.tsv
	# g1-g2-normalized_vst_counts.tsv
	# g1-g2-DEA.tsv

        ### Differencial expression analyses
        dea <- DESeq(fi, parallel = T)
        dea_contrast <- results(dea, lfcThreshold= 1, altHypothesis="greaterAbs", parallel = T)
        #save(dea_contrast, file = paste0("./DEA/by_groups/", group1,"_vs_", group2, ".Rdata", sep = ""))
        
	dea_df <- as.data.frame(dea_contrast)
	##  Save files for trinotate
	
	baseMeanA <- rowMeans(counts(dea, normalized=TRUE)[,colData(dea)$Group == group1])
        baseMeanB <- rowMeans(counts(dea, normalized=TRUE)[,colData(dea)$Group == group2])
        
        res = cbind(baseMeanA, baseMeanB, dea_df)
        res = cbind(sampleA=group1, sampleB=group2, as.data.frame(res))
	res = res[complete.cases(res),]
        res = as.data.frame(res[order(res$pvalue),])
        write.table(res, file = paste0(outdir, group1, "_vs_", group2,"/Trinity.isoform.counts.matrix.", group1, "_vs_", group2, ".DESeq.DE_results"), sep = "\t", quote = F )
	write.table(res, file = paste0(outdir,"Trinity.isoform.counts.matrix.", group1, "_vs_", group2, ".DESeq.DE_results"), sep = "\t", quote = F )

        write.table(cbind(as.character(filter$Group), filter$Sample_file), file= paste0(outdir, group1, "_vs_", group2,"/Trinity.isoform.counts.matrix.", group1, "_vs_", group2, ".DESeq.samples"), col.names = F, row.names = F, sep = "\t", quote = F)
        
        ### filter DEGs 
        filter_1 <- dea_df[complete.cases(dea_df),]
        filter_2 <- filter_1[filter_1$padj < 0.05,]
        filter_3 <- filter_2[abs(filter_2$log2FoldChange) > 1,]
	
        ## PlotMA  
        ylim <- c(-10, 10)
        png(paste0(outdir,"plotMA-", group1, "-", group2, ".png", sep=""),  width = 15*1.3, height = 15, res = 320, units = "cm", pointsize = 12, bg = "white")
        DESeq2::plotMA(dea_contrast, ylim=ylim, main = paste0( group1 ,"-", group2, sep = "" ), alpha = 0.05) 
        dev.off()
        ## volcano plot
        filter_1$test = filter_1$padj < 0.05 & abs(filter_1$log2FoldChange) > 1
        filter_1 = rownames_to_column(filter_1, var='ensgene')
        
        g = ggplot(filter_1, aes(x = log2FoldChange,
                                 y = -log10(padj), 
                                 name = ensgene)) +
                geom_point(aes(colour = test), size = 1, alpha = 0.3) +
                scale_colour_manual(values=c('black', 'red')) +
                geom_vline(xintercept = 1, colour ='green', linetype = 3) +
                geom_vline(xintercept = -1, colour ='green', linetype = 3) +
                geom_hline(yintercept = -log10(0.05), colour = 'blue', linetype = 3) +
                labs(title = paste( group1, "vs", group2)) +
                theme_bw() +
                theme( text = element_text(size = 22), 
                       panel.border = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.line = element_line(colour = "black"),
                       legend.position = 'none')
        
        ggsave( g, filename = paste0(outdir,"volcanoplot-",  group1, "-", group2, ".png"), units = "cm", width = 15*1.3, height = 15, dpi = 320)
        
	filter_3$names <- row.names(filter_3)
	write.table(filter_3, file = paste0(outdir, group1, "_vs_", group2, "/", group1, "-", group2, "-", "DEA.tsv", sep = ""), sep = "\t", quote = F, col.names = T)
	write.table(filter_3, file = paste0(outdir, group1, "-", group2, "-", "DEA.tsv", sep = ""), sep = "\t", quote = F, col.names = T)

        df_fi <- assay(fi)
        df_vst <- assay(vst)	

	write.table(df_fi, file = paste0(outdir, group1, "_vs_", group2, "/", group1, "-", group2, "-", "filtered_raw_counts.tsv", sep = ""), sep = "\t" ,quote = F, col.names = T)
	write.table(df_fi, file = paste0(outdir, group1, "-", group2, "-", "filtered_raw_counts.tsv", sep = ""), sep = "\t" ,quote = F, col.names = T)

        write.table(df_vst, file = paste0(outdir, group1, "_vs_", group2, "/", group1, "-", group2, "-", "vst_counts.tsv", sep = ""), sep = "\t" ,quote = F, col.names = T)
	write.table(df_vst, file = paste0(outdir, group1, "-", group2, "-", "vst_counts.tsv", sep = ""), sep = "\t" ,quote = F, col.names = T)

	write.table(filter, file = paste0(outdir, group1, "_vs_", group2, "/", group1, "-", group2, "-", "metadata.tsv", sep = ""), sep = "\t", quote = F, col.names = T)
	write.table(filter, file = paste0(outdir, group1, "-", group2, "-", "metadata.tsv", sep = ""), sep = "\t", quote = F, col.names = T)
}

#dea_by_group("./../data/metadata_complete.csv", "R_270_B0", "R_10_B0", "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")
#dea_by_group("./../data/metadata_complete.csv", "NR_270_B0", "NR_10_B0", "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")
#dea_by_group("./../data/metadata_complete.csv", "R_270_B", "R_10_B", "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")
#dea_by_group("./../data/metadata_complete.csv", "NR_270_B", "NR_10_B", "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")
#dea_by_group("./../data/metadata_complete.csv", "R_270_M", "R_10_M", "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")
#dea_by_group("./../data/metadata_complete.csv", "NR_270_M", "NR_10_M", "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")
#dea_by_group("./../data/metadata_complete.csv", "R_270_P", "R_10_P", "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")
#dea_by_group("./../data/metadata_complete.csv", "NR_270_P", "NR_10_P", "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")

dea_by_group("./../data/metadata_complete.csv", "R_270_B0", "NR_270_B0", "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")
dea_by_group("./../data/metadata_complete.csv", "R_10_B0", "NR_10_B0", "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")
dea_by_group("./../data/metadata_complete.csv", "R_270_B", "NR_270_B", "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")
dea_by_group("./../data/metadata_complete.csv", "R_10_B", "NR_10_B", "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")
dea_by_group("./../data/metadata_complete.csv", "R_270_M", "NR_270_M", "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")
dea_by_group("./../data/metadata_complete.csv", "R_10_M", "NR_10_M", "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")
dea_by_group("./../data/metadata_complete.csv", "R_270_P", "NR_270_P", "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")
dea_by_group("./../data/metadata_complete.csv", "R_10_P", "NR_10_P", "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")



