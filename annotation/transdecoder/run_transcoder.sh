#!/bin/bash
#$ -q all.q
#$ -cwd
#$ -V
#$ -pe smp 1

module load transdecoder/5.5.0

TransDecoder.LongOrfs -t /Storage/data1/riano/Sugarcane/NitrogenResponsiveGenotypes/Salmon//All_masked__100p__unmasked.fasta -O ./
#TransDecoder.Predict -t /Storage/data1/riano/Sugarcane/NitrogenResponsiveGenotypes/Salmon//All_masked__100p__unmasked.fasta -O ./
