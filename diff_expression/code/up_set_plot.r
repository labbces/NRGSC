# Upsetplot
#library(UpsetR)
install.packages("ggVennDiagram")
library(ggVennDiagram)
#dirs
setwd("/Storage/data1/jorge.munoz/NRGSC.new/code")
data_dir="/Storage/data1/jorge.munoz/NRGSC.new/results/gene/"
# R
# up
df1 <- read.table(paste0(data_dir, "R_270_B_vs_R_10_B_up", sep=""), col.names = "B")
df2 <- read.table(paste0(data_dir, "R_270_B0_vs_R_10_B0_up", sep=""), col.names = "B0")
df3 <- read.table(paste0(data_dir, "R_270_M_vs_R_10_M_up", sep=""), col.names = "M")
df4 <- read.table(paste0(data_dir, "R_270_P_vs_R_10_P_up", sep=""), col.names = "P")
# down
df5 <- read.table(paste0(data_dir, "R_270_B_vs_R_10_B_down", sep = ""), col.names = "B")
df6 <- read.table(paste0(data_dir, "R_270_B0_vs_R_10_B0_down", sep = ""), col.names = "B0")
df7 <- read.table(paste0(data_dir, "R_270_M_vs_R_10_M_down", sep = ""), col.names = "M")
df8 <- read.table(paste0(data_dir, "R_270_P_vs_R_10_P_down",  sep = ""), col.names = "P")
# Example lists (replace with your actual data)
up <- list(B=as.character(df1$B), B0=as.character(df2$B0), M=as.character(df3$M), P=as.character(df4$P))
down <- list(B=as.character(df5$B), B0=as.character(df6$B0), M=as.character(df7$M), P=as.character(df8$P))
# NR
# up
df9 <- read.table(paste0(data_dir, "NR_270_B_vs_NR_10_B_up", sep = ""), col.names = "B")
df10 <- read.table(paste0(data_dir, "NR_270_B0_vs_NR_10_B0_GO_up", sep= ""), col.names = "B0")
df11 <- read.table(paste0(data_dir, "NR_270_M_vs_NR_10_M_up", sep = ""), col.names = "M")
df12 <- read.table(paste0(data_dir, "NR_270_M_vs_NR_10_M_up", sep = ""), col.names = "P")
# down
df13 <- read.table(paste0(data_dir, "NR_270_B_vs_NR_10_B_down", sep = ""), col.names = "B")
df14 <- read.table(paste0(data_dir, "NR_270_B0_vs_NR_10_B0_GO_up", sep = ""), col.names = "B0")
df15 <- read.table(paste0(data_dir, "NR_270_M_vs_NR_10_M_up", sep = ""), col.names = "M")
df16 <- read.table(paste0(data_dir, "NR_270_P_vs_NR_10_P_up", sep = ""), col.names = "P")
# data
R_up <- list(B=as.character(df1$B), B0=as.character(df2$B0), M=as.character(df3$M), P=as.character(df4$P))
R_down <- list(B=as.character(df5$B), B0=as.character(df6$B0), M=as.character(df7$M), P=as.character(df8$P))
NR_up <- list(B=as.character(df9$B), B0=as.character(df10$B0), M=as.character(df11$M), P=as.character(df12$P))
NR_down <- list(B=as.character(df13$B), B0=as.character(df14$B0), M=as.character(df15$M), P=as.character(df16$P))
#upset(data_test, sets = c("set1", "set2", "set3"))
ggVennDiagram(R_up, force_upset = T)
ggVennDiagram(NR_up, force_upset = T, )
ggVennDiagram(R_down, force_upset = T)
ggVennDiagram(NR_down, force_upset = T)


