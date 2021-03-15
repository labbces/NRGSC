#!/bin/bash
#$ -q all.q
#$ -cwd
#$ -V
#$ -pe smp 4

echo $NSLOTS
module load R/4.0.0

#bash bash_scripts.sh
Rscript DEA_1v.r
Rscript venn_diagrams.r
echo termine!!!!
