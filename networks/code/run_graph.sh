#!/bin/bash
#$ -q all.q
#$ -cwd
#$ -V
#$ -pe smp 1

module load R/4.0.0

#Rscript make_expression_matrix.r
#Rscript make_graph_diff_expr.r
#Rscript make_graph_diff_corrected.r
#Rscript make_graph_from_tpm.r

#Rscript topGO.r
#Rscript make_heatmaps_for_all_modules.r
Rscript GO_DEG.r
#Rscript KEGG_DEG.r
#Rscript make_heatmaps_for_modules.r
#Rscript get_adjTOlist.r
#Rscript graph_importante2.r
