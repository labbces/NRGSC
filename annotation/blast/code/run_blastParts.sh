#!/bin/bash
#$ -t 1-100
#$ -cwd
#$ -q all.q
#$ -tc 4
#$ -pe smp 5


module load blast/2.8.1+
export m=$SGE_TASK_ID
export sprot=/Storage/data1/jorge.munoz/Trinotate-Trinotate-v3.2.2/RUN/uniprot_sprot.pep

#rm /Storage/data1/jorge.munoz/blast/blast_parts/results/*
blastx -query ./../data/parts/All_masked__100p__unmasked.part-${m}.fasta -db $sprot -num_threads 5 -max_target_seqs 6 -outfmt 6 -evalue 1e-3 >> ./../results/blastx.outfmt6.full.ok
blastp -query ./../data/parts/longest_orfs.part-${m}.pep -db $sprot -num_threads 5 -max_target_seqs 6 -outfmt 6 -evalue 1e-3 >> ./../results/blastp.outfmt6.full.ok

#ls ./../results/ | grep blastx.outfmt6* | cat > ./../results/full_blastx.outfmt6
#ls ./../results/ | grep blastp.outfmt6* | cat > ./../results/full_blastp.outfmt6

#rm ./../results/blastx.outfmt6_*
#rm ./../results/blastp.outfmt6_*

