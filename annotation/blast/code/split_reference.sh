#!/bin/bash
#$ -q all.q
#$ -cwd

ensemble="/Storage/data1/riano/Sugarcane/NitrogenResponsiveGenotypes/Salmon//All_masked__100p__unmasked.fasta"
#prot="/Storage/data1/jorge.munoz/TRANSDECODER/longest_orfs.pep"

mkdir -p ./../data/parts/

perl ./fasta-splitter.pl --n-parts 100 --out-dir ./../data/parts/ --line-length --nopad $ensemble
#perl /Storage/data1/jorge.munoz/blast/blast_parts/code/fasta-splitter.pl --n-parts 100 --out-dir ./../data/parts/ --line-length --nopad $prot
