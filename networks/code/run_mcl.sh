#!/bin/bash
#$ -t 1-12
#$ -cwd
#$ -q all.q
#$ -pe smp 1
#$ -tc 20

export m=$((SGE_TASK_ID-1))
module load mcl/14-137
LC_NUMERIC="en_US.UTF-8".
num=($(LC_NUMERIC="en_US.UTF-8". seq 1.3 0.5 6))

in=/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/pearson_abs.current.triples
out=/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/pearson_mcl/

mkdir $out
mcl $in -I ${num[${m}]} -te 1 --abc -o ${out}out.${num[${m}]} 
