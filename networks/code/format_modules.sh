#!/bin/bash
#$ -q all.q
#$ -cwd
#$ -V
#$ -pe smp 1

IN_FILE=/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/pearson_abs.current.triples
OUT_FILE=/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/pearson_mcl/out
num=($(LC_NUMERIC="en_US.UTF-8". seq 1.3 0.5 6))

## GET COUNTS IN EACH CLUSTER
for i in "${num[@]}"
do
	echo $i
	awk -F' ' '{print NF}' ${OUT_FILE}.${i} > ${OUT_FILE}.${i}.number
done

### REFORMAT MODULES OUTPUT

dir=/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/pearson_mcl
cd $dir
for file in $(ls $dir | grep -v number | grep -v formated)
do
	echo $file	
	rm -fv ${dir}/${file}.formated.csv
	count=1
	for line in $(cat $file)
	do
        	sed -n "${count}p" $file | awk -F" " -v c="$count" '{for(i=1; i<=NF; i++) {print $i,c}}' >> ${file}.formated.csv
        	let count++
	done
done
cd -
