# Upsetplot
#library(UpSetR)
#library(ggVennDiagram)
library(ggplot2)
library(wesanderson)

# set working directory
setwd("~/upset/")
# R
# up
df1 <- read.table("R_270_B_vs_R_10_B_up", col.names = "B")
df2 <- read.table("R_270_B0_vs_R_10_B0_up" , col.names = "B0")
df3 <- read.table("R_270_M_vs_R_10_M_up", col.names = "M")
df4 <- read.table("R_270_P_vs_R_10_P_up",  col.names = "P")
# down
df5 <- read.table("R_270_B_vs_R_10_B_down", col.names = "B")
df6 <- read.table("R_270_B0_vs_R_10_B0_down" , col.names = "B0")
df7 <- read.table("R_270_M_vs_R_10_M_down", col.names = "M")
df8 <- read.table("R_270_P_vs_R_10_P_down",  col.names = "P")
# Example lists 
up <- list(B=as.character(df1$B), B0=as.character(df2$B0), M=as.character(df3$M), P=as.character(df4$P))
down <- list(B=as.character(df5$B), B0=as.character(df6$B0), M=as.character(df7$M), P=as.character(df8$P))
# NR
# up
df9 <- read.table("NR_270_B_vs_NR_10_B_up", col.names = "B")
df10 <- read.table("NR_270_B0_vs_NR_10_B0_GO_up" , col.names = "B0")
df11 <- read.table("NR_270_M_vs_NR_10_M_up", col.names = "M")
df12 <- read.table("NR_270_M_vs_NR_10_M_up",  col.names = "P")
# down
df13 <- read.table("NR_270_B_vs_NR_10_B_down", col.names = "B")
df14 <- read.table("NR_270_B0_vs_NR_10_B0_GO_down" , col.names = "B0")
df15 <- read.table("NR_270_M_vs_NR_10_M_down", col.names = "M")
df16 <- read.table("NR_270_P_vs_NR_10_P_down",  col.names = "P")
# data
R_up <- list(B=as.character(df1$B), B0=as.character(df2$B0), M=as.character(df3$M), P=as.character(df4$P))
R_down <- list(B=as.character(df5$B), B0=as.character(df6$B0), M=as.character(df7$M), P=as.character(df8$P))
NR_up <- list(B=as.character(df9$B), B0=as.character(df10$B0), M=as.character(df11$M), P=as.character(df12$P))
NR_down <- list(B=as.character(df13$B), B0=as.character(df14$B0), M=as.character(df15$M), P=as.character(df16$P))





# all 
# R
R <- list(up_B1=as.character(df1$B), up_B0=as.character(df2$B0), up_M=as.character(df3$M), up_P=as.character(df4$P), down_B1=as.character(df5$B), down_B0=as.character(df6$B0), down_M=as.character(df7$M), down_P=as.character(df8$P))
# NR
NR <- list(up_B1=as.character(df9$B), up_B0=as.character(df10$B0), up_M=as.character(df11$M), up_P=as.character(df12$P), down_B1=as.character(df13$B), down_B0=as.character(df14$B0), down_M=as.character(df15$M), down_P=as.character(df16$P))

# all all
all <- list(R_up_B1=as.character(df1$B), R_up_B0=as.character(df2$B0), R_up_M=as.character(df3$M), R_up_P=as.character(df4$P), R_down_B1=as.character(df5$B), R_down_B0=as.character(df6$B0), R_down_M=as.character(df7$M), R_down_P=as.character(df8$P), NR_up_B1=as.character(df9$B), NR_up_B0=as.character(df10$B0), NR_up_M=as.character(df11$M), NR_up_P=as.character(df12$P), NR_down_B1=as.character(df13$B), NR_down_B0=as.character(df14$B0), NR_down_M=as.character(df15$M), NR_down_P=as.character(df16$P))

# Function to convert list to binary matrix for UpSetR
list_to_matrix <- function(input_list) {
  # Get all unique elements across all sets
  all_elements <- unique(unlist(input_list))
  
  # Create binary matrix
  matrix_data <- matrix(0, nrow = length(all_elements), ncol = length(input_list))
  rownames(matrix_data) <- all_elements
  colnames(matrix_data) <- names(input_list)
  
  # Fill matrix with 1s where element is present in set
  for (i in seq_along(input_list)) {
    matrix_data[input_list[[i]], i] <- 1
  }
  
  return(matrix_data)
}

# Convert your existing lists to matrices
R_up_matrix <- list_to_matrix(R_up)
R_down_matrix <- list_to_matrix(R_down)
NR_up_matrix <- list_to_matrix(NR_up)
NR_down_matrix <- list_to_matrix(NR_down)
R_matrix <- list_to_matrix(R)
NR_matrix <- list_to_matrix(NR)
all_matrix <- list_to_matrix(all)

# Now you can use UpSetR's upset function
library(UpSetR)

# Create upset plots
R_up_upset <- upset(as.data.frame(R_up_matrix), 
                    sets = colnames(R_up_matrix),
                    order.by = "freq")

R_down_upset <- upset(as.data.frame(R_down_matrix), 
                      sets = colnames(R_down_matrix),
                      order.by = "freq")

NR_up_upset <- upset(as.data.frame(NR_up_matrix), 
                     sets = colnames(NR_up_matrix),
                     order.by = "freq")

NR_down_upset <- upset(as.data.frame(NR_down_matrix), 
                       sets = colnames(NR_down_matrix),
                       order.by = "freq")

# For larger matrices, you might want to limit intersections shown
R_upset <- upset(as.data.frame(R_matrix), 
                 sets = colnames(R_matrix),
                 nsets = ncol(R_matrix),
                 nintersects = 30,
                 order.by = "freq")

NR_upset <- upset(as.data.frame(NR_matrix), 
                  sets = colnames(NR_matrix),
                  nsets = ncol(NR_matrix),
                  nintersects = 30,
                  order.by = "freq")

all_upset <- upset(as.data.frame(all_matrix), 
                   sets = colnames(all_matrix),
                   nsets = ncol(all_matrix),
                   nintersects = 40,
                   order.by = "freq")


colors <- wes_palette("Darjeeling1")
bar_colors <- rep("black", 103)  # Default color for all bars
bar_colors[14] <- colors[1]  # Highlight bar 14 in red

sets_ordered <- rev(c("NR_down_B0", "NR_down_B1", "NR_down_M", "NR_down_P",
                  "NR_up_B0", "NR_up_B1", "NR_up_M", "NR_up_P",
                  "R_down_B0", "R_down_B1", "R_down_M", "R_down_P",
                  "R_up_B0", "R_up_B1", "R_up_M", "R_up_P"))


png("simple_pretty_upset_plot.png", width = 1600*2.2, height = 1000*2.2, res = 300)
upset(as.data.frame(all_matrix), 
      sets = sets_ordered,
      keep.order = T,
      nsets = ncol(all_matrix),
      nintersects = 20,
      order.by = "freq",
      decreasing = TRUE,
      mb.ratio = c(0.6, 0.4),
      text.scale = c(1.5, 1.2, 1.2, 1, 1.5, 2.5),
      point.size = 3,
      line.size = 1.5,
      mainbar.y.label = paste0("Shared Differentially Expressed" , "\n", "Groups of Transcripts"),
      sets.x.label = paste0("Quantity of Differentialy Expressed Groups" ,  "\n", "of Transcripts in Each Contrast"),
      sets.bar.color = "black",
      main.bar.color = bar_colors,  # Use custom color vector
      matrix.color = "black",
      shade.color = "black")
dev.off()



######################### UPSET FOR CORRELATION ###############################

# read data
clorophyll_a <- read.table("clorophyll_a_bonferroni_005.txt", col.names = "Clorophyll_a")
clorophyll_b <- read.table("clorophyll_b_bonferroni_005.txt", col.names = "Clorophyll_b")
pepcase <- read.table("pepcase_bonferroni_005.txt", col.names = "Pepcase")
rubisco <- read.table("rubisco_bonferroni_005.txt", col.names = "Rubisco")
genotype <- read.table("genotype_bonferroni_005.txt", col.names = "Genotype")
total_clorophyll <- read.table("total_clorophyll_bonferroni_005.txt", col.names = "Total_Clorophyll")
condition <- read.table("condtition_bonferroni_005.txt", col.names = "Condition")
# make a list with all the data frames
all_cor <- list(clorophyll_a = as.character(clorophyll_a$clorophyll_a), clorophyll_b = as.character(clorophyll_b$clorophyll_b), pepcase = as.character(pepcase$pepcase), rubisco = as.character(rubisco$rubisco), genotype = as.character(genotype$Genotype), total_clorophyll = as.character(total_clorophyll$Total_Clorophyll), condition = as.character(condition$Condition))
# transform list to binary matrix

cor_matrix <- list_to_matrix(all_cor)

png("simple_pretty_cor_upset_plot.png", width = 1600*1.2, height = 1000*1.5, res = 300)
upset(as.data.frame(cor_matrix), 
      keep.order = T,
      nsets = ncol(cor_matrix),
      nintersects = 20,
      order.by = "freq",
      decreasing = TRUE,
      mb.ratio = c(0.7, 0.3),
      text.scale = c(1.2, 1.2, 1, 1, 1.2, 1.2),
      point.size = 3,
      line.size = 1.5,
      mainbar.y.label = paste0("Shared Correlated Modules"),
      sets.x.label = paste0("Quantity of Modules",  "\n", "Correlated to Features"),
      sets.bar.color = "black",
      main.bar.color = "black",  
      matrix.color = "black",
      shade.color = "black")
dev.off()
