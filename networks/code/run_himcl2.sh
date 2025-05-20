#!/bin/bash
#$ -q all.q
#$ -cwd
#$ -V
#$ -pe smp 1

module load R
module load hipmcl
module load miniconda3
conda activate himcl


IN_FILE=/Storage/data1/jorge.munoz/NRGSC.new/networks/results/abs.current.triples
OUT_FILE=/Storage/data1/jorge.munoz/NRGSC.new/networks/results/out.hipmlc
#export OMP_NUM_THREADS=1


## EVALUATE INFLATION VALUE
#rm -f /Storage/data1/jorge.munoz/NRGSC.new/networks/results/I_test.csv
#LC_ALL=C seq 0.1 0.1 7 > param
#for i in $(cat param)
#do
#mpirun -np $NSLOTS hipmcl -M $IN_FILE -I ${i} -per-process-mem 30 -o $OUT_FILE${i}
#No=$(echo $OUT_FILE_I${i} | wc -l)
#echo ${No},${i} >> /Storage/data1/jorge.munoz/NRGSC.new/networks/results/I_test.csv
#done
#conda deactivate
#Rscript inflation_plot.r

## RUN INFLATIN SELECTIONED
mpirun -np 1 hipmcl -M $IN_FILE -I 2 -per-process-mem 100 -o $OUT_FILE
### GET COUNTS IN EACH CLUSTER
awk -F' ' '{print NF}' $OUT_FILE > number.${OUT_FILE}
### REFORMAT MODULES OUTPUT
rm -f /Storage/data1/jorge.munoz/NRGSC.new/networks/results/modules_formated_ok.csv
count=1
for line in $(cat $OUT_FILE)
do
  	sed -n "${count}p" $OUT_FILE | awk -F" " -v c="$count" '{for(i=1; i<=NF; i++) {print $i,c}}' >> /Storage/data1/jorge.munoz/NRGSC.new/networks/results/modules_formated_ok.csv
        let count++
done

