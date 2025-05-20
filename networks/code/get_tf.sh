# SCRIPT TO GET TRANSCRIPTION FACTORS IN THE CO-EXPRESSION NETWORK
# Define variables 
TFs=/Storage/data1/jorge.munoz/GET_TFS/results/TF_no_Orphans.ids
NET_IDs=/Storage/data1/jorge.munoz/NRGSC.new/networks/data/network.id
OUT=/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/TF/new/TFs_in_network.ids
# TF list is per peptide, we need TF listo per "gene"
cut -d"." -f1 /Storage/data1/jorge.munoz/GET_TFS/results/TF_no_Orphans.ids | uniq | awk '{print $0".gen"}' > /Storage/data1/jorge.munoz/GET_TFS/results/TF_no_Orphans_gene.ids
TFs_ok=/Storage/data1/jorge.munoz/GET_TFS/results/TF_no_Orphans_gene.ids

# make some comm magic
comm -12 <(sort $TFs_ok) <(sort $NET_IDs) > $OUT
#grep -w -f $NET_IDs $TFs_ok > $OUT
# GET hub genes
module load R/4.0.0
Rscript get_hub.r
# GET TF in hub genes
HUB=/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/TF/new/newhub_genes.txt
TF_NET=/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/TF/new/TFs_in_network.ids
TF_HUB=/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/TF/new/TFs_hubs.id
comm -12 <(sort $HUB) <(sort $TF_NET) > $TF_HUB
# Now we search which peptides belong to the transcripts
TF_TMP=/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/TF/new/TFs_hubs_no_gene.id
cut -f1 -d"." $TF_HUB > $TF_TMP
# see wich peptides are TF
TF_OUT=/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/TF/new/TFs_out.id
TF_FULL=/Storage/data1/jorge.munoz/GET_TFS/results/TF_no_Orphans.TXT
grep -f $TF_TMP $TF_FULL > $TF_OUT



