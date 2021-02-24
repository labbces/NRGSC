main_folder="./"
mkdir -p libraries vst rld RUVs RUVs/R RUVs/NR RUVg RUVg/R RUVg/NR RUVr RUVr/R RUVr/NR
## get samples names in the folder
ls $main_folder | grep Sample_ >  "$main_folder"/samples2.txt
## build txgen.txt
#cut -f1 "$main_folder"/Sample_19/quant.sf > tmp
#paste tmp tmp > "$main_folder"/txgen2.txt
#rm tmp
# remove first line
#sed -i -e "1d"  "$main_folder"/txgen2.txt
# write first line
#sed -i '1s/^/transcript_id	gene_id\n/' "$main_folder"/txgen2.txt
##

