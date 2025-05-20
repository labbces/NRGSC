# Read nitrogen modules
Nmods <- read.table("/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/pearson_GO_1.8/all_mods.txt", header = F)
colnames(Nmods) <- "Mod No"

# files
#modules_path <- "./../../results/modules_formated.csv"
modules_path <- "/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/pearson_mcl/out.1.8.formated.csv"
vst_path <- "/Storage/data1/jorge.munoz/NRGSC.new/networks/data/all_vst_counts.tsv"
metadata_path <-"/Storage/data1/jorge.munoz/NRGSC.new/networks/data/metadata_complete.csv"

#libraries
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

modules <- read.table(modules_path, row.names = 1, header = F)
colnames(modules) <- c("module_No")
# read full vst matrix 
#vst <- read.table(vst_path)
metadata <- read.table(metadata_path, sep = ",", header = T)

trait_genotype <-as.numeric(c("0", "0", "0", "0", "0", "0", "0", "0", "1", "1", "1", "1", "1", "1", "1", "1"))
trait_nitrogen_condition <-as.numeric(c("1", "1", "1", "1","0", "0", "0", "0","1", "1", "1", "1","0", "0", "0", "0"))
trait_leaf_section <- as.numeric(c("1", "2", "3", "4", "1", "2", "3", "4", "1", "2", "3", "4", "1", "2", "3", "4"))
trait_R_up_response <- as.numeric(c("0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "1", "1", "1", "1"))
trait_interaction_NG <- trait_genotype * trait_nitrogen_condition
# 10_NR_B0 10_NR_B 10_NR_M 10_NR_P 270_NR_B0 270_NR_B 270_NR_M 270_NR_P 10_R_B0 10_R_B 10_R_M 10_R_P 270_R_BO 270_R_B 270_R_M 270_R_P 

trait_total_clorophyll <- as.numeric(c("0.34", "0.58", "0.92", "0.92", "0.82", "2.45", "2.84", "1.59", "0.43", "0.97", "1.19", "1.32", "1.04", "2.46", "3.05", "2.43"))
trait_clorophyll_a <- as.numeric(c("0.126", "0.324", "0.504", "0.522", "0.45", "1.323", "1.395", "0.81", "0.18", "0.522", "0.675", "0.792", "0.648", "1.395", "1.764", "1.08"))
trait_clorophyll_b<- as.numeric(c("0.2136", "0.2616", "0.35", "0.40", "0.3816", "1.03", "1.29", "0.79", "0.2568", "0.335", "0.52", "0.51", "0.3864", "1.1", "1.25", "1.43"))
trait_rubisco <- as.numeric(c("2", "4", "4.72", "5.47", "2.1", "8.22", "9.49", "7.01", "1.41", "2.45", "3.11", "4.72", "1.36", "9.92", "7.75", "4.29"))
trait_pepcase <- as.numeric(c("0.62", "1.37", "1.37", "2.55", "2.47", "6.36", "9.58", "4.27", "0.39", "0.91", "1", "0.84", "0.97", "3.58", "8.58", "4.81"))



# 1 = B0
# 2 = B
# 3 = M
# 4 = P
# get_cor_module(module10_dat, trait_genotype, "M10")
get_cor_module <- function(x,trait, module){
df <- x
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
pca <- prcomp(t(df2))
eigen <-pca$x[,1]
cor <- cor.test(eigen, trait, method = "spearman", exact=FALSE)
pvalue <- cor$p.value
rho <- cor$estimate

df_out <- data.frame(rho=rho,
p = pvalue,
Module = module,
row.names = NULL
)
return(df_out)
}

condition_table <- data.frame(matrix(ncol = 3, nrow = dim(Nmods)[1]))
genotype_table <- data.frame(matrix(ncol = 3, nrow = dim(Nmods)[1]))
leaf_segment_table <- data.frame(matrix(ncol = 3, nrow = dim(Nmods)[1]))
up_R_response_table <- data.frame(matrix(ncol = 3, nrow = dim(Nmods)[1]))

total_clorophyl_table <- data.frame(matrix(ncol = 3, nrow = dim(Nmods)[1]))
clorophyll_a_table <- data.frame(matrix(ncol = 3, nrow = dim(Nmods)[1]))
clorophyll_b_table <- data.frame(matrix(ncol = 3, nrow = dim(Nmods)[1]))
rubisco_table <- data.frame(matrix(ncol = 3, nrow = dim(Nmods)[1]))
pepcase_table <- data.frame(matrix(ncol = 3, nrow = dim(Nmods)[1]))
interaction_NG_table <- data.frame(matrix(ncol = 3, nrow = dim(Nmods)[1]))

colnames(condition_table) <- c("rho", "pvalue", "module")
colnames(genotype_table) <- c("rho", "pvalue", "module")
colnames(leaf_segment_table) <- c("rho", "pvalue", "module")
colnames(up_R_response_table) <- c("rho", "pvalue", "module")

colnames(total_clorophyl_table) <- c("rho", "pvalue", "module")
colnames(clorophyll_a_table) <- c("rho", "pvalue", "module")
colnames(clorophyll_b_table) <- c("rho", "pvalue", "module")
colnames(rubisco_table) <- c("rho", "pvalue", "module")
colnames(pepcase_table) <- c("rho", "pvalue", "module")
colnames(interaction_NG_table) <- c("rho", "pvalue", "module")

for (i in Nmods$'Mod No'){
condition_table[i,] <- get_cor_module(eval(as.name(paste0("module",i,"_dat"))), trait_nitrogen_condition, paste0("M", i))
genotype_table[i,] <- get_cor_module(eval(as.name(paste0("module",i,"_dat"))), trait_genotype, paste0("M", i))
leaf_segment_table [i,] <- get_cor_module(eval(as.name(paste0("module",i,"_dat"))), trait_leaf_section, paste0("M", i))
up_R_response_table[i,] <- get_cor_module(eval(as.name(paste0("module",i,"_dat"))), trait_R_up_response, paste0("M", i))

total_clorophyl_table[i,] <- get_cor_module(eval(as.name(paste0("module",i,"_dat"))), trait_total_clorophyll, paste0("M", i))
clorophyll_a_table[i,] <- get_cor_module(eval(as.name(paste0("module",i,"_dat"))), trait_clorophyll_a, paste0("M", i))
clorophyll_b_table[i,] <- get_cor_module(eval(as.name(paste0("module",i,"_dat"))), trait_clorophyll_b, paste0("M", i))
rubisco_table[i,] <- get_cor_module(eval(as.name(paste0("module",i,"_dat"))),trait_rubisco , paste0("M", i))
pepcase_table[i,] <- get_cor_module(eval(as.name(paste0("module",i,"_dat"))),trait_pepcase , paste0("M", i))
interaction_NG_table[i,] <- get_cor_module(eval(as.name(paste0("module",i,"_dat"))), trait_interaction_NG, paste0("M", i))


}
## Print full matrices
# Create a list of all trait tables
trait_tables <- list(
		       condition = condition_table,
		         genotype = genotype_table,
			   leaf_segment = leaf_segment_table,
			     R_response = up_R_response_table,
			       total_chlorophyll = total_clorophyl_table,
			         chlorophyll_a = clorophyll_a_table,
				   chlorophyll_b = clorophyll_b_table,
				     rubisco = rubisco_table,
				       pepcase = pepcase_table,
				         interaction_NG = interaction_NG_table
		       
		     )

# Initialize empty matrices
n_modules <- nrow(Nmods)
n_traits <- length(trait_tables)
correlation_matrix <- matrix(nrow = n_modules, ncol = n_traits)
pvalue_matrix <- matrix(nrow = n_modules, ncol = n_traits)

# Set row and column names
rownames(correlation_matrix) <- paste0("M", 1:n_modules)
rownames(pvalue_matrix) <- paste0("M", 1:n_modules)
colnames(correlation_matrix) <- names(trait_tables)
colnames(pvalue_matrix) <- names(trait_tables)

# Fill matrices
for (i in 1:length(trait_tables)) {
	  trait_name <- names(trait_tables)[i]
  correlation_matrix[,i] <- trait_tables[[i]]$rho
    pvalue_matrix[,i] <- trait_tables[[i]]$pvalue
  
}

# Round values for better readability
correlation_matrix <- round(correlation_matrix, 6)
pvalue_matrix <- round(pvalue_matrix, 6)

# Save matrices to files (optional)
write.csv(correlation_matrix, "correlation_matrix.csv")
write.csv(pvalue_matrix, "pvalue_matrix.csv")

# Print preview of matrices
print("Correlation Matrix:")
print(head(correlation_matrix))
print("\nP-value Matrix:")
print(head(pvalue_matrix))


##
filtered_condition <- mutate(condition_table, adjusted = p.adjust(condition_table$pvalue, method = "bonferroni")) %>% mutate(condition_table, fdr = p.adjust(condition_table$pvalue, method = "fdr")) %>% filter(adjusted  < 0.05 & abs(rho) > 0.7)
filtered_genotype <- mutate(genotype_table, adjusted = p.adjust(genotype_table$pvalue, method = "bonferroni")) %>% mutate(genotype_table, fdr = p.adjust(genotype_table$pvalue, method = "fdr")) %>% filter(adjusted  < 0.05 & abs(rho) > 0.7)
filtered_leaf_segment <- mutate(leaf_segment_table, adjusted = p.adjust(leaf_segment_table$pvalue, method = "bonferroni")) %>% mutate(genotype_table, fdr = p.adjust(genotype_table$pvalue, method = "fdr")) %>% filter(adjusted  < 0.05 & abs(rho) > 0.5)
filtered_interaction_NG <- mutate(interaction_NG_table, adjusted = p.adjust(interaction_NG_table$pvalue, method = "bonferroni")) %>% mutate(interaction_NG_table, fdr = p.adjust(interaction_NG_table$pvalue, method = "fdr")) %>% filter(adjusted < 0.05 & abs(rho) > 0.7)

filtered_pepcase <- mutate(pepcase_table, adjusted = p.adjust(pepcase_table$pvalue, method = "bonferroni")) %>% mutate(pepcase_table, fdr = p.adjust(pepcase_table$pvalue, method = "fdr")) %>% filter(adjusted < 0.05 & abs(rho) > 0.7)
filtered_rubisco <- mutate(rubisco_table, adjusted = p.adjust(rubisco_table$pvalue, method = "bonferroni")) %>% mutate(rubisco_table, fdr = p.adjust(rubisco_table$pvalue, method = "fdr")) %>% filter(adjusted < 0.05 & abs(rho) > 0.7)
filtered_total_clorophyl <- mutate(total_clorophyl_table, adjusted = p.adjust(total_clorophyl_table$pvalue, method = "bonferroni")) %>% mutate(total_clorophyl_table, fdr = p.adjust(total_clorophyl_table$pvalue, method = "fdr")) %>% filter(adjusted < 0.05 & abs(rho) > 0.7)

filtered_clorophyl_a <- mutate(clorophyll_a_table, adjusted = p.adjust(clorophyll_a_table$pvalue, method = "bonferroni")) %>% mutate(clorophyll_a_table, fdr = p.adjust(clorophyll_a_table$pvalue, method = "fdr")) %>% filter(adjusted < 0.05 & abs(rho) > 0.7)
filtered_clorophyl_b <- mutate(clorophyll_b_table, adjusted = p.adjust(clorophyll_b_table$pvalue, method = "bonferroni")) %>% mutate(clorophyll_b_table, fdr = p.adjust(clorophyll_b_table$pvalue, method = "fdr")) %>% filter(adjusted < 0.05 & abs(rho) > 0.7)

#filtered <- filter(table, pvalue < 0.01 & abs(rho) > 0.7)
#write.table(filtered_condition, "../../results/perNCONDITION/module_nitrogen_condition_correlation_fdr_0.5.csv", row.names = F, quote = F) 
#write.table(filtered_genotype, "../../results/perNCONDITION/module_genotype_correlation_fdr_0.5.csv", row.names = F, quote = F)
#write.table(filtered_leaf_segment, "../../results/perNCONDITION/module_leaf_segment_fdr_0.5.csv", row.names = F, quote = F)
# Write tables
write.table(filtered_condition$module, "../../results/perNCONDITION/correlation/condtition_bonferroni_005.txt", row.names = F, quote = F, col.names = F)
write.table(filtered_genotype$module, "../../results/perNCONDITION/correlation/genotype_bonferroni_005.txt", row.names = F, quote = F, col.names = F)
write.table(filtered_pepcase$module, "../../results/perNCONDITION/correlation/pepcase_bonferroni_005.txt", row.names = F, quote = F, col.names = F)
write.table(filtered_rubisco$module, "../../results/perNCONDITION/correlation/rubisco_bonferroni_005.txt", row.names = F, quote = F, col.names = F)
write.table(filtered_total_clorophyl$module, "../../results/perNCONDITION/correlation/total_clorophyll_bonferroni_005.txt", row.names = F, quote = F, col.names = F)
write.table(filtered_clorophyl_a$module, "../../results/perNCONDITION/correlation/clorophyll_a_bonferroni_005.txt", row.names = F, quote = F, col.names = F)
write.table(filtered_clorophyl_b$module, "../../results/perNCONDITION/correlation/clorophyll_b_bonferroni_005.txt", row.names = F, quote = F, col.names = F)

## fix to run in server
# activate R-env in conda and run as is
setwd("/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/correlation")
library(ggVennDiagram)
library(ggplot2)
df1 <- read.table("total_clorophyll_bonferroni_005.txt", col.names = "Module")
df2 <- read.table("clorophyll_a_bonferroni_005.txt", col.names = "Module")
df3 <- read.table("clorophyll_b_bonferroni_005.txt", col.names = "Module")
df4 <- read.table("condtition_bonferroni_005.txt", col.names = "Module")
df5 <- read.table("genotype_bonferroni_005.txt", col.names = "Module")
df6 <- read.table("pepcase_bonferroni_005.txt", col.names = "Module")
df7 <- read.table("rubisco_bonferroni_005.txt", col.names = "Module")

all <- list(total_chlorophyll = df1$Module, chlorophyll_a = df2$Module, chlorophyll_b = df3$Module, nitrogen_condition = df4$Module, genotype = df5$Module, pepcase = df6$Module, rubisco = df7$Module)
all_upset_plot <- ggVennDiagram(all, force_upset = T ,order.set.by = "none", nintersects = 17)
ggsave("upsetplot_correlations.png", all_upset_plot, "png", dpi = 300)
