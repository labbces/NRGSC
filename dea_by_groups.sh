#!/bin/bash
#$ -q all.q
#$ -cwd
#$ -V
#$ -pe smp 10
echo $NSLOTS
module load R/4.0.0
Rscript DEA_1v_cluster_by_groups.r
echo termine!!!!
