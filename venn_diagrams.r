install.packages("VennDiagram", 
lib = "./libraries",
repos = 'http://cran.us.r-project.org')

# Load library
library(VennDiagram, lib.loc =  "./libraries")

# Generate  sets
setwd("./DEA")
setB_NR <- read.table(file = "filter_3NR_270_B_vs_NR_10_B.csv", sep = ",", header = T)
setB0_NR <- read.table(file = "filter_3NR_270_B0_vs_NR_10_B0.csv", sep = ",", header = T)
setB_R <- read.table(file = "filter_3R_270_B_vs_R_10_B.csv", sep = ",", header = T)
setB0_R <- read.table(file = "filter_3R_270_B0_vs_R_10_B0.csv", sep = ",", header = T)
setM_R <- read.table(file = "filter_3R_270_M_vs_R_10_M.csv", sep = ",", header = T)
setM_NR <- read.table(file = "filter_3NR_270_M_vs_NR_10_M.csv", sep = ",", header = T)
setP_NR <- read.table(file = "filter_3NR_270_P_vs_NR_10_P.csv", sep = ",", header = T)
setP_R <- read.table(file = "filter_3R_270_P_vs_R_10_P.csv", sep = ",", header = T)
# get column of with names of genes deferentially expressed
set_NR_B <- setB_NR$names
set_NR_B0 <- setB0_NR$names
set_R_B <- setB_R$names
set_R_B0 <- setB0_R$names
set_R_M <- setM_R$names
set_NR_M <- setM_NR$names
set_NR_P <- setP_NR$names
set_R_P <- setP_R$names
# make venn diagrams for interesting comparations , max 4 groups
# venn 1
venn.diagram(filename = "test.venn1", x = list(set_R_B, set_R_B0, set_R_M, set_R_P),
             category.names = c("set_R_P", "set_NR_P", "set_NR_M", "set_R_M"), height = 6000, width = 6000)
# venn 2
venn.diagram(filename = "test.venn2", x = list(set_NR_B, set_NR_B0, set_NR_M, set_R_P),
             category.names = c("set_R_P", "set_NR_P", "set_NR_M", "set_R_M"), height = 6000, width = 6000)
# venn 3
venn.diagram(filename = "test.venn3", x = list(set_NR_B, set_R_B0, set_NR_P, set_R_M),
             category.names = c("set_R_P", "set_NR_P", "set_NR_M", "set_R_M"), height = 6000, width = 6000)
